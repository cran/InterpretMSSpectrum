#'@title findMAIN.
#'
#'@aliases findMAIN
#'
#'@description
#'\code{findMAIN} will evaluate an ESI spectrum for the potential main adducts, rank obtained suggestions and allow the deduction of the neutral mass of the measured molecule.
#'
#'@details
#' Electrospray ionization (ESI) mass spectra frequently contain a number of different adduct ions, multimers and in-source fragments ([M+H]+, [M+Na]+, [2M+H]+, [M+H-H2O]+), making it difficult to decide on the compound's neutral mass.
#' This functions aims at determining the main adduct ion and its type (protonated, sodiated etc.) of a spectrum, allowing subsequent database searches e.g. using MS-FINDER, SIRIUS or similar.
#'
#' @param spec A mass spectrum. Either a matrix or data frame, the first two columns of which are assumed to contain the 'mz' and 'intensity' values, respectively.
#' @param adductmz Manually specified peak for which \code{adducthyp} should be tested, or 'NULL' (default), to test all main peaks. What is a main peak, is governed by \code{mainpkthr}.
#' @param ionmode Ionization mode, either "positive" or "negative". Can be abbreviated.
#' @param adducthyp Adduct hypotheses to test for each main peak. Defaults to \code{c("[M+H]+","[M+Na]+","[M+K]+")} for positive mode and \code{c("[M-H]-","[M+Cl]-","[M+HCOOH-H]-")}.
#' @param ms2spec Second spectrum limiting main peak selection. If available, MS^E or bbCID spectra may allow further exclusion of false positive adduct ions, as ions of the intact molecule (protonated molecule, adduct ions) should have lower intensity in the high-energy trace than in low-energy trace.
#' @param rules Adduct/fragment relationships to test, e.g. \code{c("[M+Na]+", "[M+H-H2O]+")}, or 'NULL' for default set (see \code{\link{Adducts}})
#' @param mzabs Allowed mass error, absolute (Da).
#' @param ppm Allowed mass error, relative (ppm), which is _added_ to 'mzabs'.
#' @param mainpkthr Intensity threshold for main peak selection, relative to base peak.
#' @param collapseResults If a neutral mass hypothesis was found more than once (due to multiple adducts suggesting the same neutral mass), return only the one with the highest adduct peak. Should normally kept at \code{TRUE}, the default.
#'
#' @return A list-like 'findMAIN' object for which 'print', 'summary' and 'plot' methods are available.
#' 
#' @references Jaeger C, Meret M, Schmitt CA, Lisec J (2017), RCM, accepted.
#'
#' @examples
#' \donttest{
#' utils::data(esi_spectrum, package = "InterpretMSSpectrum")
#' fmr <- findMAIN(esi_spectrum)
#' plot(fmr)
#' head(summary(fmr))
#' InterpretMSSpectrum(fmr[[1]], precursor=263, ionization="ESI+", formula_db="ESI.db")
#'}
#'
#' @export

findMAIN <-
  function(spec,
           adductmz = NULL,
           ionmode = c("positive", "negative")[1],
           adducthyp = NULL,
           ms2spec = NULL,
           rules = NULL,
           mzabs = 0.01,
           ppm = 5,
           mainpkthr = 0.005,
           collapseResults = TRUE) {
    Debug = FALSE
    getlabel <- function(x, mhyp) {
      x
    }
    nummatch <- function(x, table, interval, min.only = TRUE) {
      d <- abs(x - table)
      d[d > interval] <- NA
      out <- if (min.only)
        which.min(d)
      else
        which(!is.na(d))
      return(if (length(out) > 0)
        out
        else
          NA)
    }
    checkSpec <- function(spec) {
      if (!inherits(spec, c("matrix", "data.frame")))
        stop("invalid spectrum")
      s <- spec
      s <- s[order(s[, 1]), , drop = FALSE]
      s <- s[!duplicated(s[, 1]), , drop = FALSE]
      s[, 2][is.na(s[, 2]) | s[, 2] < 0] <- 0
      ## s[,2] <- s[,2] / max(s[,2]) * 100
      if (any(is.na(match(
        c("isogr", "iso", "charge"), colnames(s)
      )))) {
        ## get isotope annotation if not available
        s <-
          findiso(s,
                  mzabs = mzabs,
                  intthr = 0.03,
                  CAMERAlike = TRUE)
      }
      return(s)
    }
    scaleSpec <- function(spec) {
      s <- spec
      s[, 2] <- s[, 2] / max(s[, 2]) * 100
      return(s)
    }
    getMainPeaks <- function(spec, intthr, ms2spec = NULL) {
      s <- spec
      if (is.null(ms2spec)) {
        good <- (s[, "iso"] == 0 | is.na(s[, "iso"])) &
          s[, 2] > (max(s[, 2]) * intthr)
      } else {
        ms2int <- ms2spec[, 2][sapply(
          ms2spec[, 1],
          nummatch,
          table = s[, 1],
          interval = 0.01,
          min.only = TRUE
        )]
        good <- (s[, "iso"] == 0 | is.na(s[, "iso"])) &
          s[, 2] > (max(s[, 2]) * intthr) & s[, 2] > ms2int
      }
      mainpks <-
        s[, 1][good][order(s[, 2][good], decreasing = TRUE)]
      return(mainpks)
    }
    getDefaultRules <- function(ionmode) {
      Adducts <- NULL
      utils::data(Adducts, envir = environment(), package = "InterpretMSSpectrum") # Adducts
      on.exit(rm(Adducts))
      rules <- switch(ionmode, positive = Adducts$Positive, negative = Adducts$Negative)
      if (is.null(rules)) stop("unknown ionmode")
      return(rules)
    }
    generateRules <- function(rules) {
      getRuleFromIonSymbol(rules)
    }
    predictPeaksFromRuleset <- function(neutral_mass, ruleset) {
      r <- ruleset
      return((neutral_mass * r[, "nmol"] + r[, "massdiff"]) / abs(r[, "charge"]))
    }
    resolveConflicts <- function(s, ruleset, rules.found) {
      allPeaksUnique <- all(sapply(rules.found, length) > 1)
      allRulesUnique <-
        all(!duplicated(rules.found) & !is.na(rules.found))
      if (allPeaksUnique &&
          allRulesUnique)
        return(rules.found)
      ## first resolve double peak assignments (same rule assigned to more than one peak)
      if (!allPeaksUnique) {
        notok <- which(sapply(rules.found, length) > 1)
        for (i in notok) {
          conflictingPks.Intensities <- s[rules.found[[i]], 2]
          rules.found[[i]] <-
            rules.found[[i]][which.max(conflictingPks.Intensities)]
        }
      }
      ## then resolve double rule assignments (same peak assigned multiple rules (mostly charge-related)
      if (any(idx <-
              duplicated(rules.found) & !is.na(rules.found))) {
        notok <- which(idx)
        nonuniquePks <-
          unique(unlist(rules.found[which(idx)]))
        for (pk.idx in nonuniquePks) {
          ##pk.idx <- rules.found[[i]]
          r.idx <- which(!is.na(unlist(rules.found)) &
                           unlist(rules.found) == pk.idx)
          rch <- ruleset[, "charge"][r.idx]
          pkch <- s[, "charge"][pk.idx]
          ##message("r.idx:", r.idx, "; rch:", rch, "; pk.idx:", pk.idx, "; pkch:", pkch)
          r.idx.wrong <-
            ifelse(!is.na(pkch), r.idx[which(rch != pkch)],
                   r.idx[which(rch != which.min(rch))])
          for (k in 1:length(r.idx.wrong))
            rules.found[[r.idx.wrong[k]]] <- NA
        }
      }
      return(rules.found)
    }
    findMatchingRules <- function(s, prec, adducthyp = NULL, ruleset, mzabs, ppm) {
        rules.found <- vector("numeric", nrow(s))
        neutral_mass <-
          if (is.null(adducthyp))
            prec
        else
          prec - adducthyp
        expectedPeaks <-
          predictPeaksFromRuleset(neutral_mass, ruleset)
        test.idx <- which(is.na(s[, "iso"]) | s[, "iso"] == 0)
        prec.idx <- NULL
        if (!is.null(adducthyp)) {
          prec.idx <- nummatch(prec, s[, 1], mzabs + ppm * prec / 1e6, min.only =
                                 T)
          prec.idx.extended <-
            nummatch(prec, s[, 1], mzabs + ppm * prec / 1e6, min.only = F)
          if (is.numeric(prec.idx.extended))
            test.idx <- test.idx [!test.idx %in% prec.idx.extended]
        }
        test.mz <- vector("numeric", length = nrow(s))
        test.mz[test.idx] <- s[, 1][test.idx]
        rules.found <- lapply(expectedPeaks, function(x) {
          nummatch(x,
                   table = test.mz,
                   mzabs + ppm * x / 1e6,
                   min.only = F)
        })
        rules.found <- resolveConflicts(s, ruleset, rules.found)
        rules.found <- unlist(rules.found) # now unique
        r.idx <- which(!is.na(rules.found))
        ##exclude.r.idx <- r.idx[nummatch(adducthyp, ruleset[,4][r.idx], 0.001, F)]
        ##if(is.numeric(exclude.r.idx)) r.idx <- r.idx[-exclude.r.idx]
        pk.idx <- rules.found[r.idx]
        dmz <- abs(s[pk.idx, , drop = F][, 1] - expectedPeaks[r.idx])
        return(list(r.idx, pk.idx, dmz, prec.idx))
      }
    scoreMatchingRules <- function(s,
                                   ruleset,
                                   matchingRules,
                                   maxExplainedAdducts,
                                   maxExplainedIntensity) {
      out <- matrix(NA, ncol = 7, nrow = 0)
      colnames(out) <- c(
        "adducts_explained",
        "medppm",
        "int_perc",
        "mass_score",
        "int_score",
        "supp_isos",
        "total_score"
      )
      if (Debug)
        out <-
        cbind(out, matrix(
          ncol = 4,
          nrow = 0,
          dimnames = list(NULL, c(
            "score1", "score2", "score3", "score4"
          ))
        ))
      if (is.null(s) ||
          (length(matchingRules[[1]]) + length(matchingRules[[4]])) == 0)
        return(rbind(out, NA))
      e1 = 0.5
      e2 = 0.2
      w1 = 0.5
      w2 = 0
      w3 = 0.1
      w4 = 0.4
      
      minppm = 2
      isoTrue = 1
      isoFalse = 0
      isoNeutral = 0.75
      s.deiso <- s[is.na(s[, 4]) | s[, 4] == 0,]
      pk.idx <- c(matchingRules[[2]], matchingRules[[4]])
      r.idx <- matchingRules[[1]]
      dmz <- c(matchingRules[[3]], if (is.null(matchingRules[[4]])) NULL else NA)
      adducts_explained <- length(c(matchingRules[[2]], matchingRules[[4]]))
      ## ppm score
      ppmVals <- dmz / s[pk.idx, , drop = F][, 1] * 1e6
      medppm <- stats::median(ppmVals, na.rm = TRUE)
      ppmVals[ppmVals < minppm] <- minppm
      ppmVals <- (minppm / ppmVals) ^ e2
      ppmVals1 <- ppmVals
      ppmVals1[is.na(ppmVals1)] <- (minppm / max(minppm * 2, ppm / 2)) ^ e2
      ## charge score
      isoCharge <- s[, 5][pk.idx]
      ruleCharge <- c(ruleset[r.idx, ][, 3], if (is.null(matchingRules[[4]])) NULL else 1)
      isNA <- is.na(isoCharge)
      isoChargeVals <- rep(isoNeutral, length(pk.idx))
      isoChargeVals[!isNA & isoCharge == ruleCharge] <- 1
      isoChargeVals[!isNA & isoCharge != ruleCharge] <- 0
      supp_isos <- length(which(isoCharge == ruleCharge))
      ## intensity score
      intVals <- s[, 2][pk.idx] / sum(s.deiso[, 2]) / maxExplainedIntensity
      intExpl <- sum(s[, 2][pk.idx] / sum(s.deiso[, 2]))
      ## total score
      score1 <- sum(((intVals ^ e1 / sum(intVals ^ e1)) * intExpl) * isoChargeVals)
      score2 <- sum(((intVals ^ e1 / sum(intVals ^ e1)) * intExpl) * ppmVals1)
      score3 <- sum(isoChargeVals) / adducts_explained
      score4 <- adducts_explained / maxExplainedAdducts
      total_score <- score1 * w1 + score2 * w2 + score3 * w3 + score4 * w4
      res <- cbind(
        adducts_explained = adducts_explained,
        medppm = medppm,
        int_perc = intExpl,
        mass_score = score3,
        int_score = score1,
        supp_isos = supp_isos,
        total_score = total_score
      )
      if (Debug)
        res <- cbind(
          res,
          score1 = score1,
          score2 = score2,
          score3 = score3,
          score4 = score4
        )
      out <- rbind(out, res)
      return(out)
    }
    formatResults <- function(s, matchingRules, ruleset, adductmz, adducthyp, scores) {
        r.idx <- matchingRules[[1]]
        pk.idx <- c(matchingRules[[2]], matchingRules[[4]])
        deltamz <- c(matchingRules[[3]], NA)
        adductname <- c(ruleset[r.idx, "name"],
                        ifelse(is.null(names(adducthyp)),
                               round(adducthyp, 4),
                               names(adducthyp)))
        ppm <- abs(deltamz) / s[, 1][pk.idx] * 1e6
        s.out <- cbind(s,
                       adduct = NA,
                       ppm = NA,
                       label = NA)
        s.out[pk.idx, ][, "adduct"] <- adductname
        s.out[pk.idx, ][, "ppm"] <- ppm
        s.out[, "label"] <- s.out[, "adduct"]
        scores.out <-
          cbind(
            adductmz = adductmz,
            adducthyp = adducthyp,
            neutral_mass = adductmz - adducthyp,
            scores
          )
        attr(s.out, "scores") <- scores.out
        s.out
      }
    getExplainedIntensity <- function(s, matchingRules) {
      s.deiso <- s[is.na(s[, 4]) | s[, 4] == 0, , drop = F]
      isAdduct <- c(matchingRules[[2]], matchingRules[[4]])
      sum(s[, 2][isAdduct] / sum(s.deiso[, 2]))
    }
    collapseResultSet <- function(ResultSet) {
      if (length(ResultSet) < 2)
        return(ResultSet)
      rs <- ResultSet
      x <-
        data.frame(matrix(ncol = ncol(attr(rs[[1]], "scores")),
                          nrow = length(rs)))
      colnames(x) <- colnames(attr(rs[[1]], "scores"))
      for (i in 1:length(rs))
        x[i,] <- attr(rs[[i]], "scores")
      x[, "idx"] <- 1:nrow(x)
      mztab <- rs[[1]][, 1]
      inttab <- rs[[1]][, 2]
      x[, "int"] <-
        inttab[sapply(x[, "adductmz"], function(mz)
          nummatch(mz, mztab, 0.0001))]
      x[, "nm_grp"] <-
    stats::cutree(stats::hclust(stats::dist(x[, "neutral_mass"])), h = 0.015)
      cs <- tapply(x[, "total_score"], x[, "nm_grp"], max)
      x[, "cs"] <- cs[match(x[, "nm_grp"], names(cs))]
      x <- x[rev(order(x[, "cs"], x[, "int"])), ]
      x <- x[!duplicated(x[, "nm_grp"]), ]
      x[, "total_score"] <- x[, "cs"]
      for (i in 1:nrow(x))
        attr(rs[[x[i, "idx"]]], "scores") <-
        x[i, 1:ncol(attr(rs[[1]], "scores"))]
      idx.good <- x[, "idx"]
      return(rs[idx.good])
    }
    ##
    ## main
    ##
    s <- s.in <- spec
    if (nrow(s) == 0 || sum(s[, 2]) == 0) {
      warning("spectrum without peaks")
      return(NULL)
    }
    s <- checkSpec(s)
    if (!is.null(ms2spec))
      ms2spec <- checkSpec(ms2spec)
    ionmode <- match.arg(ionmode, c("positive", "negative"))
    if (is.null(adducthyp))
      adducthyp <- switch(
        ionmode,
        positive = c("[M+H]+", "[M+Na]+", "[M+K]+"),
        negative = c("[M-H]-", "[M+Cl]-", "[M+HCOOH-H]-")
      )
    adducthyp <-
      unlist(sapply(adducthyp, getRuleFromIonSymbol)[4, ])
    if (nrow(s) == 1) {
      ## return something useful for spectra with only one line
      warning(sprintf(
        "single-line spectrum - wild guessing %s",
        round(adducthyp[1]),
        4
      ))
      s.out <- data.frame(
        s[, 1:5, drop = FALSE],
        adduct = NA,
        mDa = NA,
        ppm = NA,
        label = NA
      )
      score <- scoreMatchingRules(NULL)
      attr(s.out, "scores") <-
        data.frame(
          adductmz = s[, 1],
          adducthyp = adducthyp[1],
          neutral_mass = s[, 1] - adducthyp[1],
          score
        )
      out <- list(s.out)
      class(out) <- "findMAIN"
      attr(out, "adducthyp_tested") <- names(adducthyp)[1]
      attr(out, "rules_tested") <- NULL
      return(out)
    }
    if (nrow(s) == 2)
      adducthyp <- adducthyp[1]
    # test only first adducthyp for poor specs
    if (is.null(rules))
      rules <- getDefaultRules(ionmode)
    rules <- generateRules(rules)
    prec <- adductmz
    if (is.null(prec)) {
      prec <- getMainPeaks(s, intthr = mainpkthr, ms2spec = ms2spec)
    }
    prectab <- expand.grid(prec = prec, adducthyp = adducthyp)
    prectab <-
      cbind(prectab, neutral_mass = prectab[, 1] - prectab[, 2])
    s <- scaleSpec(s)
    if (!is.null(ms2spec))
      ms2spec <- scaleSpec(ms2spec)
    matchingRules <- lapply(1:nrow(prectab), function(i) {
      findMatchingRules(s,
                        prectab[i, 1],
                        prectab[i, 2],
                        rules,
                        mzabs = mzabs,
                        ppm = ppm)
    })
    maxExplainedAdducts <-
      max(sapply(lapply(matchingRules, "[[", 1), length) + 1)
    maxExplainedIntensity <-
      max(sapply(matchingRules, getExplainedIntensity, s = s), na.rm = T)
    scores <- lapply(1:length(matchingRules), function(i)
      scoreMatchingRules(
        s,
        rules,
        matchingRules[[i]],
        maxExplainedAdducts,
        maxExplainedIntensity
      ))
    out <- lapply(1:length(matchingRules), function(i)
      formatResults(s, matchingRules[[i]], rules,
                    prectab[i, ][, 1], prectab[i, ][, 2], scores[[i]]))
    if (collapseResults)
      out <- collapseResultSet(out)
    scores <-
      round(sapply(out, function(x)
        attr(x, "scores")[, "total_score"]), 2)
    adducthyps <-
      sapply(out, function(x)
        attr(x, "scores")[, "adducthyp"])
    adductrank <-
      sapply(adducthyps, nummatch, table = adducthyp, interval = 0.001)
    # order by specified adduct order
    o <- order(scores,-adductrank, decreasing = TRUE)
    out <- out[o]
    attr(out, "adducthyp_tested") <- names(adducthyp)
    attr(out, "rules_tested") <- rules[, 1]
    class(out) <- "findMAIN"
    return(out)
  }


print.findMAIN <-
  function (x)
  {
    nres <- length(x)
    i = 1
    scores <- summary(x)
    scores1 <- attr(x[[i]], "scores")
    nprec <- length(unique(sapply(x, function(x)
      attr(x, "scores")[,
                        "adductmz"])))
    nadduct <-
      length(unique(sapply(x, function(x)
        attr(x, "scores")[,
                          "adducthyp"])))
    message(
      sprintf(
        "Analyzed %d neutral mass hypotheses (%d peaks * %d adducts), kept %d",
        nprec * nadduct,
        nprec,
        nadduct,
        nres
      )
    )
    message(
      sprintf(
        "Selected %.4f as [M%+.2f], neutral %.4f, score %.2f",
        scores1[, 1],
        scores1[, 2],
        scores1[, 3],
        scores1[, 10]
      )
    )
    print(x[[i]])
  }


plot.findMAIN <-
  function (x,
            rank = 1,
            correct_mass = NULL,
            ...)
  {
    if (length(rank) > 1) {
      opar <- graphics::par(mfrow = grDevices::n2mfrow(length(rank)))
      graphics::par(mar = c(2, 2, 2, 1))
      on.exit(graphics::par(opar))
    }
    idx <- rank
    if (any(idx > length(x))) {stop("object shorter than requested numer of ranks")}
    lidx <- length(idx)
    legend_text_col <- matrix(1, nrow = ncol(attr(x[[1]], "scores")), ncol = lidx)
    if (lidx > 1) {
      sm <- summary(x)[idx, ]
      legend_text_col[7, which(sm$mass_score == max(sm$mass_score))] <- 3
      legend_text_col[8, which(sm$int_score == max(sm$int_score))] <- 3
      legend_text_col[9, which(sm$supp_isos == max(sm$supp_isos))] <- 3
    }
    for (i in 1:lidx) {
      cols <- x[[idx[i]]][,"charge"]
      cols[is.na(cols)] <- 0
      cols <- cols+1
      cols[cols>=3] <- 6
      cols[!is.na(x[[idx[i]]][,"label"])] <- 4
      InterpretMSSpectrum::PlotSpec(x = x[[idx[i]]], cutoff = 0, cols = cols, txt = x[[idx[i]]][,c("mz","label")], ...)
      graphics::legend("topright", legend = paste(colnames(attr(x[[idx[i]]], "scores")), round(as.numeric(attr(x[[idx[i]]], "scores")), 2)), bty = "n", cex = 0.75, text.col = legend_text_col[, i])
      mhyp <- attr(x[[idx[i]]], "scores")[, "neutral_mass"]
      color <- if (!is.null(correct_mass)) {
        if (abs(mhyp - correct_mass) < 0.01) 3 else 2
      } else {
        1
      }
      mtext(sprintf("[%d] %.4f", idx[i], round(mhyp, 4)), col = color)
    }
  }

summary.findMAIN <-
  function (x)
  {
    attrname <- "scores"
    out <- data.frame(matrix(ncol = ncol(attr(x[[1]], attrname)), nrow = length(x)))
    colnames(out) <- colnames(attr(x[[1]], attrname))
    for (i in 1:length(x)) {
      out[i,] <- attr(x[[i]], attrname)
    }
    for (i in c(1, 2, 3)) out[, i] <- round(out[, i], 4)
    for (i in c(5, 6, 7, 8, 10))  out[, i] <- round(out[, i], 2)
    return(out)
  }
