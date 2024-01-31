#' @title findMAIN.
#'
#' @aliases findMAIN
#'
#' @description \code{findMAIN} will evaluate an ESI spectrum for the potential main adducts, 
#'     rank obtained suggestions and allow the deduction of the neutral mass of the measured 
#'     molecule.
#'
#' @details Electrospray ionization (ESI) mass spectra frequently contain a number of different 
#'     adduct ions, multimers and in-source fragments \code{[M+H]+, [M+Na]+, [2M+H]+, [M+H-H2O]+}, 
#'     making it difficult to decide on the compound's neutral mass. This functions aims 
#'     at determining the main adduct ion and its type (protonated, sodiated etc.) of a spectrum, 
#'     allowing subsequent database searches e.g. using MS-FINDER, SIRIUS or similar.     
#'
#' @param spec A mass spectrum. Either a matrix or data frame, the first two columns of which 
#'     are assumed to contain the 'mz' and 'intensity' values, respectively.
#' @param adductmz Manually specified peak for which \code{adducthyp} should be tested, 
#'     or 'NULL' (default), to test all main peaks. What is a main peak, is governed by 
#'     \code{mainpkthr}.
#' @param ionmode Ionization mode, either "positive" or "negative". Can be abbreviated.
#' @param adducthyp Adduct hypotheses to test for each main peak. Defaults to 
#'     \code{c("[M+H]+","[M+Na]+","[M+K]+")} for positive mode and 
#'     \code{c("[M-H]-","[M+Cl]-","[M+HCOOH-H]-")} for negative mode.
#' @param ms2spec Second spectrum limiting main peak selection. If available, MS^E or bbCID 
#'     spectra may allow further exclusion of false positive adduct ions, as ions of the intact 
#'     molecule (protonated molecule, adduct ions) should have lower intensity in the high-energy 
#'     trace than in low-energy trace.
#' @param rules Adduct/fragment relationships to test, e.g. \code{c("[M+Na]+", "[M+H-H2O]+")}, 
#'     or 'NULL' for default set (see \code{\link{Adducts}})
#' @param mzabs Allowed mass error, absolute (Da).
#' @param ppm Allowed mass error, relative (ppm), which is _added_ to 'mzabs'.
#' @param mainpkthr Intensity threshold for main peak selection, relative to base peak.
#' @param collapseResults If a neutral mass hypothesis was found more than once (due to multiple 
#'     adducts suggesting the same neutral mass), return only the one with the highest adduct peak. 
#'     Should normally kept at \code{TRUE}, the default.
#'
#' @return A list-like 'findMAIN' object for which 'print', 'summary' and 'plot' methods are available.
#'     Each list element represents a potential spectra annotation, ranked according to 
#'     a combined score. The spectrum is annotated with columns indicating the determined
#'     isotopic groups (isogr) and their likely charge. Further, information on the potential 
#'     set of adducts and their ppm error is attached.
#'     The score aims to integrate all this information using formula S=sum(w_i x s_i).
#'     In short, we sum up i weighted score components (currently i=4). Currently these
#'     components are calculated based on the explained intensity (adduct sets which
#'     annotate a higher amount of the total spectrum intensity are better), the mass error
#'     (adduct sets with lower mass error are better), the support by isotopic peaks
#'     (adduct sets with fitting isotopes are better) and the number of adducts (adduct 
#'     sets with a larger number of adducts are better).
#'     The individual scores for each adduct set are attached as an attribute to the
#'     respective list element and can be easily observed by applying the 'summary' or 
#'     the 'plot' function on the 'findMAIN' object.
#' 
#' @references Jaeger C, Meret M, Schmitt CA, Lisec J (2017), <doi:10.1002/rcm.7905>.
#'
#' @examples
#' \donttest{
#' utils::data(esi_spectrum, package = "InterpretMSSpectrum")
#' fmr <- InterpretMSSpectrum::findMAIN(esi_spectrum)
#' plot(fmr)
#' head(summary(fmr))
#' InterpretMSSpectrum::InterpretMSSpectrum(fmr[[1]], precursor=263, param="ESIpos")
#' fmr <- InterpretMSSpectrum::findMAIN(esi_spectrum[6:9,], adducthyp = "[M+H]+")
#' plot(fmr)
#' 
#' # set up a spectrum containing a double charged peak
#' spec <- data.frame(mz = c(372.1894, 372.6907, 373.1931, 380), int = c(100, 40, 8, 2))
#' InterpretMSSpectrum:::findiso(spec)
#' # allow a double charged adduct hypothesis (not standard)
#' fmr <- InterpretMSSpectrum::findMAIN(spec, adducthyp = c("[M+H]+", "[M+2H]2+"))
#' summary(fmr)
#' attr(fmr[[1]],"scores")
#' plot(fmr, rank = 1:4)
#' plot(fmr, rank = 2)
#' 
#' # add the correct M+H to this spectrum as a minor peak
#' spec <- rbind(spec, c(742.3648+1.007, 10))
#' (fmr <- InterpretMSSpectrum::findMAIN(spec, adducthyp = c("[M+H]+", "[M+2H]2+")))
#' summary(fmr)
#' plot(fmr, rank = 1)
#' plot(fmr, rank = 2)
#' 
#' # compare specific hypotheses manually
#' # get correct result
#' InterpretMSSpectrum::findMAIN(spec, adductmz = 743.3718, adducthyp = "[M+H]+")
#' # enforce wrong result
#' InterpretMSSpectrum::findMAIN(spec, adductmz = 743.3718, adducthyp = "[M+2H]2+")
#' }
#'
#' @export

findMAIN <- function(
  spec,
  adductmz = NULL,
  ionmode = c("positive", "negative")[1],
  adducthyp = NULL,
  ms2spec = NULL,
  rules = NULL,
  mzabs = 0.01,
  ppm = 5,
  mainpkthr = 0.005,
  collapseResults = TRUE
) {
  Debug <- FALSE
  getlabel <- function(x, mhyp) {
    x
  }
  nummatch <- function(x, table, interval, min.only = TRUE) {
    d <- abs(x - table)
    d[d > interval] <- NA
    out <- if (min.only) {
      which.min(d)
    } else {
      which(!is.na(d))
    }
    return(if (length(out) > 0) {
      out
    } else {
      NA
    })
  }
  checkSpec <- function(spec) {
    if (!inherits(spec, c("matrix", "data.frame"))) {
      stop("invalid spectrum")
    }
    s <- spec
    s <- s[order(s[, 1]), , drop = FALSE]
    s <- s[!duplicated(s[, 1]), , drop = FALSE]
    s[, 2][is.na(s[, 2]) | s[, 2] < 0] <- 0
    ## s[,2] <- s[,2] / max(s[,2]) * 100
    if (any(is.na(match(c("isogr", "iso", "charge"), colnames(s))))) {
      ## get isotope annotation if not available
      s <- findiso(s, mzabs = mzabs, intthr = 0.03, CAMERAlike = TRUE)
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
    rules <- switch(ionmode,
      positive = Adducts$Positive,
      negative = Adducts$Negative
    )
    if (is.null(rules)) stop("unknown ionmode")
    rules <- getRuleFromIonSymbol(rules)
    return(rules)
  }
  predictPeaksFromRuleset <- function(neutral_mass, ruleset) {
    r <- ruleset
    return((neutral_mass * r[, "nmol"] + r[, "massdiff"]) / abs(r[, "charge"]))
  }
  resolveConflicts <- function(s, ruleset, rules.found) {
    allPeaksUnique <- all(sapply(rules.found, length) > 1)
    allRulesUnique <- all(!duplicated(rules.found) & !is.na(rules.found))
    if (allPeaksUnique && allRulesUnique) {
      return(rules.found)
    }
    ## first resolve double peak assignments (same rule assigned to more than one peak)
    if (!allPeaksUnique) {
      notok <- which(sapply(rules.found, length) > 1)
      for (i in notok) {
        conflictingPks.Intensities <- s[rules.found[[i]], 2]
        rules.found[[i]] <- rules.found[[i]][which.max(conflictingPks.Intensities)]
      }
    }
    ## then resolve double rule assignments (same peak assigned multiple rules (mostly charge-related)
    if (any(idx <- duplicated(rules.found) & !is.na(rules.found))) {
      notok <- which(idx)
      nonuniquePks <- unique(unlist(rules.found[which(idx)]))
      for (pk.idx in nonuniquePks) {
        ## pk.idx <- rules.found[[i]]
        r.idx <- which(!is.na(unlist(rules.found)) &
          unlist(rules.found) == pk.idx)
        rch <- ruleset[, "charge"][r.idx]
        pkch <- s[, "charge"][pk.idx]
        r.idx.wrong <- if (!is.na(pkch) && length(unique(rch)) > 1) { # in case peak charge is known, exclude rules of different charge
          r.idx[which(rch != pkch)]
        } else if (length(unique(rch)) > 1) { # in case peak charge is unknown but rules differ in charge exclude higher charges
          r.idx[which(rch > min(rch))]
        } else { # all rule charges are equal so just exclude all but the first
          r.idx[2:length(r.idx)]
        }
        for (k in 1:length(r.idx.wrong)) rules.found[[r.idx.wrong[k]]] <- NA
      }
    }
    return(rules.found)
  }
  findMatchingRules <- function(s, mz_test, dmz_adduct = NULL, rules, mzabs, ppm) {
    rules.found <- vector("numeric", nrow(s))
    # [JL]  on 20240126 neutral_mass calculation was put into a dedicated function, making the below modification obsolete
    # # [JL] to account for "[Mx]2+" hypotheses we need to adjust the neutral mass calculation for this double charge to be correct in finding expectedPeaks 
    # # [JL] 2 lines modified on 20230622; set fac=1 to revert this change
    # fac <- ifelse(names(dmz_adduct) %in% rules[,"name"], rules[rules[,"name"]==names(dmz_adduct),"charge"], 1)
    # neutral_mass <- ifelse(is.null(dmz_adduct), mz_test, fac * mz_test - dmz_adduct)
    # browser()
    neutral_mass <- getNeutralMass(m = mz_test, adduct_name = names(dmz_adduct), adduct_rules = rules)
    expectedPeaks <- predictPeaksFromRuleset(neutral_mass, rules)
    # cbind(rules, expectedPeaks)
    test.idx <- which(is.na(s[, "iso"]) | s[, "iso"] == 0)
    prec.idx <- NULL
    if (!is.null(dmz_adduct)) {
      prec.idx <- nummatch(mz_test, s[, 1], mzabs + ppm * mz_test / 1e6, min.only = T)
      prec.idx.extended <- nummatch(mz_test, s[, 1], mzabs + ppm * mz_test / 1e6, min.only = F)
      if (is.numeric(prec.idx.extended)) {
        test.idx <- test.idx[!test.idx %in% prec.idx.extended]
      }
    }
    test.mz <- vector("numeric", length = nrow(s))
    test.mz[test.idx] <- s[, 1][test.idx]
    rules.found <- lapply(expectedPeaks, function(x) {
      nummatch(x, table = test.mz, mzabs + ppm * x / 1e6, min.only = F)
    })
    rules.found <- resolveConflicts(s, rules, rules.found)
    rules.found <- unlist(rules.found) # now unique
    r.idx <- which(!is.na(rules.found))
    ## exclude.r.idx <- r.idx[nummatch(dmz_adduct, rules[,4][r.idx], 0.001, F)]
    ## if(is.numeric(exclude.r.idx)) r.idx <- r.idx[-exclude.r.idx]
    pk.idx <- rules.found[r.idx]
    dmz <- abs(s[pk.idx, , drop = F][, 1] - expectedPeaks[r.idx])
    return(list(r.idx, pk.idx, dmz, prec.idx))
  }
  scoreMatchingRules <- function(s, ruleset, matchingRules, maxExplainedAdducts, maxExplainedIntensity, adduct_hyp_charge = 1) {
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
    if (Debug) {
      out <-
        cbind(out, matrix(
          ncol = 4,
          nrow = 0,
          dimnames = list(NULL, c(
            "score1", "score2", "score3", "score4"
          ))
        ))
    }
    if (is.null(s) ||
      (length(matchingRules[[1]]) + length(matchingRules[[4]])) == 0) {
      return(rbind(out, NA))
    }
    e1 <- 0.5
    e2 <- 0.2
    w1 <- 0.5
    w2 <- 0
    w3 <- 0.1
    w4 <- 0.4

    minppm <- 2
    isoTrue <- 1
    isoFalse <- 0
    isoNeutral <- 0.75
    s.deiso <- s[is.na(s[, 4]) | s[, 4] == 0, ]
    pk.idx <- c(matchingRules[[2]], matchingRules[[4]])
    r.idx <- matchingRules[[1]]
    # [JL] 2 lines modified on 20230622
    #dmz <- c(matchingRules[[3]], if (is.null(matchingRules[[4]])) NULL else NA)
    dmz <- matchingRules[[3]]
    if (!is.null(matchingRules[[4]])) dmz <- c(dmz, NA)
    adducts_explained <- length(c(matchingRules[[2]], matchingRules[[4]]))
    ## ppm score
    ppmVals <- dmz / s[pk.idx, , drop = F][, 1] * 1e6
    medppm <- stats::median(ppmVals, na.rm = TRUE)
    ppmVals[ppmVals < minppm] <- minppm
    ppmVals <- (minppm / ppmVals)^e2
    ppmVals1 <- ppmVals
    ppmVals1[is.na(ppmVals1)] <- (minppm / max(minppm * 2, ppm / 2))^e2
    ## charge score
    isoCharge <- s[, 5][pk.idx]
    # [JL] 2 lines modified on 20230622
    #ruleCharge <- c(ruleset[r.idx, ][, 3], if (is.null(matchingRules[[4]])) NULL else 1)
    ruleCharge <- abs(ruleset[r.idx, ][, 3])
    # [JL] for a double charged molecule as prec we need to set the default value respectively
    if (!is.null(matchingRules[[4]])) ruleCharge <- c(ruleCharge, abs(adduct_hyp_charge))
    isNA <- is.na(isoCharge)
    isoChargeVals <- rep(isoNeutral, length(pk.idx))
    isoChargeVals[!isNA & isoCharge == ruleCharge] <- 1
    isoChargeVals[!isNA & isoCharge != ruleCharge] <- 0
    supp_isos <- length(which(isoCharge == ruleCharge))
    ## intensity score
    intVals <- s[, 2][pk.idx] / sum(s.deiso[, 2]) / maxExplainedIntensity
    intExpl <- sum(s[, 2][pk.idx] / sum(s.deiso[, 2]))
    ## total score
    score1 <- sum(((intVals^e1 / sum(intVals^e1)) * intExpl) * isoChargeVals)
    score2 <- sum(((intVals^e1 / sum(intVals^e1)) * intExpl) * ppmVals1)
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
    if (Debug) {
      res <- cbind(
        res,
        score1 = score1,
        score2 = score2,
        score3 = score3,
        score4 = score4
      )
    }
    out <- rbind(out, res)
    return(out)
  }
  formatResults <- function(s, matchingRules, ruleset, adductmz, adducthyp, scores) {
    r.idx <- matchingRules[[1]]
    pk.idx <- c(matchingRules[[2]], matchingRules[[4]])
    deltamz <- c(matchingRules[[3]], NA)
    adductname <- c(
      ruleset[r.idx, "name"],
      ifelse(is.null(names(adducthyp)),
        round(adducthyp, 4),
        names(adducthyp)
      )
    )
    ppm <- abs(deltamz) / s[, 1][pk.idx] * 1e6
    s.out <- cbind(s,
      adduct = NA,
      ppm = NA,
      label = NA
    )
    s.out[pk.idx, ][, "adduct"] <- adductname
    s.out[pk.idx, ][, "ppm"] <- ppm
    s.out[, "label"] <- s.out[, "adduct"]
    #fac <- ifelse(names(adducthyp) %in% rules[,"name"], rules[rules[,"name"]==names(adducthyp),"charge"], 1)
    #browser()
    scores.out <- cbind(
      adductmz = adductmz,
      adducthyp = adducthyp,
      #neutral_mass = fac * adductmz - adducthyp,
      neutral_mass = getNeutralMass(m = adductmz, adduct_name = names(adducthyp), adduct_rules = ruleset),
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
    if (length(ResultSet) < 2) { return(ResultSet) }
    rs <- ResultSet
    x <- data.frame(matrix(
      ncol = ncol(attr(rs[[1]], "scores")),
      nrow = length(rs),
      dimnames=list(NULL, colnames(attr(rs[[1]], "scores")))
    ))
    for (i in 1:length(rs)) {
      x[i, ] <- attr(rs[[i]], "scores")
      #x[i, "adducthyp"] <- rownames(attr(rs[[i]], "scores"))
    }
    x[, "idx"] <- 1:nrow(x)
    mztab <- rs[[1]][, 1]
    inttab <- rs[[1]][, 2]
    x[, "int"] <- inttab[sapply(x[, "adductmz"], function(mz) {
      nummatch(mz, mztab, 0.0001)
    })]
    x[, "nm_grp"] <- stats::cutree(stats::hclust(stats::dist(x[, "neutral_mass"])), h = 0.015)
    cs <- tapply(x[, "total_score"], x[, "nm_grp"], max)
    x[, "cs"] <- cs[match(x[, "nm_grp"], names(cs))]
    x <- x[rev(order(x[, "cs"], x[, "int"])), ]
    x <- x[!duplicated(x[, "nm_grp"]), ]
    x[, "total_score"] <- x[, "cs"]
    for (i in 1:nrow(x)) {
      keep_rn <- rownames(attr(rs[[x[i, "idx"]]], "scores"))
      attr(rs[[x[i, "idx"]]], "scores") <- x[i, 1:ncol(attr(rs[[1]], "scores"))]
      rownames(attr(rs[[x[i, "idx"]]], "scores")) <- keep_rn
    }
    idx.good <- x[, "idx"]
    return(rs[idx.good])
  }
  ##
  ## main
  ##
  s <- spec
  if (nrow(s) == 0 || sum(s[, 2]) == 0) {
    warning("spectrum without peaks")
    return(NULL)
  }
  s <- checkSpec(s)
  if (!is.null(ms2spec)) {
    ms2spec <- checkSpec(ms2spec)
  }
  ionmode <- match.arg(ionmode, c("positive", "negative"))
  if (is.null(adducthyp)) {
    adducthyp <- switch(ionmode,
      positive = c("[M+H]+", "[M+Na]+", "[M+K]+"),
      negative = c("[M-H]-", "[M+Cl]-", "[M+HCOOH-H]-")
    )
  }
  # [JL] ensure that 'adducthyp' becomes a named vector independent on the number of hypoteses (failed for n=1)
  adducthyp <- getRuleFromIonSymbol(adducthyp)
  adducthyp <- stats::setNames(object = adducthyp[,"massdiff"], nm = adducthyp[,"name"])
  # establish rules
  if (is.null(rules)) {
    rules <- getDefaultRules(ionmode)
  } else {
    rules <- getRuleFromIonSymbol(rules)
  }
  if (nrow(s) == 1) {
    ## return something useful for spectra with only one line
    warning(sprintf("single-line spectrum - wild guessing %s", round(adducthyp[1]), 4))
    s.out <- data.frame(
      s[, 1:5, drop = FALSE],
      adduct = NA,
      mDa = NA,
      ppm = NA,
      label = NA
    )
    score <- scoreMatchingRules(NULL)
    #fac <- unname(sapply(names(adducthyp), function(x) { generateRules(x)[,"charge"] }))
    attr(s.out, "scores") <- data.frame(
      adductmz = s[, 1],
      adducthyp = names(adducthyp)[1],
      #neutral_mass = fac[1] * s[, 1] - adducthyp[1],
      neutral_mass = getNeutralMass(m = s[, 1], adduct_name = names(adducthyp)[1], adduct_rules = rules),
      score
    )
    out <- list(s.out)
    class(out) <- "findMAIN"
    attr(out, "adducthyp_tested") <- names(adducthyp)[1]
    attr(out, "rules_tested") <- NULL
    return(out)
  }
  if (nrow(s) == 2) {
    ## test only first adducthyp for poor specs
    adducthyp <- adducthyp[1]
  }
  prec <- adductmz
  if (is.null(prec)) {
    prec <- getMainPeaks(s, intthr = mainpkthr, ms2spec = ms2spec)
  }
  prectab <- expand.grid(prec = prec, adducthyp = adducthyp)
  prectab <- cbind(prectab, neutral_mass = prectab[, 1] - prectab[, 2])
  s <- scaleSpec(s)
  if (!is.null(ms2spec)) {
    ms2spec <- scaleSpec(ms2spec)
  }
  matchingRules <- lapply(1:nrow(prectab), function(i) {
    findMatchingRules(
      s = s, 
      mz_test = prectab[i, 1], 
      dmz_adduct = prectab[i, 2], 
      rules = rules,
      mzabs = mzabs,
      ppm = ppm
    )
  })
  maxExplainedAdducts <- max(sapply(lapply(matchingRules, "[[", 1), length) + 1)
  maxExplainedIntensity <- max(sapply(matchingRules, getExplainedIntensity, s = s), na.rm = T)
  scores <- lapply(1:length(matchingRules), function(i) {
    scoreMatchingRules(
      s,
      rules,
      matchingRules[[i]],
      maxExplainedAdducts,
      maxExplainedIntensity,
      adduct_hyp_charge = rules[rules["name"]==names(prectab[i, 2]),"charge"]
    )
  })
  #browser()
  out <- lapply(1:length(matchingRules), function(i) {
    formatResults(s, matchingRules[[i]], rules, prectab[i, ][, 1], prectab[i, ][, 2], scores[[i]])
  })
  if (collapseResults) {
    out <- collapseResultSet(out)
  }
  scores <- round(sapply(out, function(x) {
    attr(x, "scores")[, "total_score"]
  }), 2)
  adducthyps <- sapply(out, function(x) {
    attr(x, "scores")[, "adducthyp"]
  })
  adductrank <- sapply(adducthyps, nummatch, table = adducthyp, interval = 0.001)
  # order by specified adduct order
  o <- order(scores, -adductrank, decreasing = TRUE)
  out <- out[o]
  attr(out, "adducthyp_tested") <- names(adducthyp)
  attr(out, "rules_tested") <- rules[, 1]
  class(out) <- "findMAIN"
  return(out)
}

#'@rdname findMAIN
#'@param x Object of class findMAIN.
#'@param rank Rank of the suggestion to plot (can be a numeric vector).
#'@param correct_mass If provided will indicate correct suggestion by green color.
#'@param ... Further plotting parameters.
#'@export
#'@method plot findMAIN
plot.findMAIN <- function (x, rank = 1, correct_mass = NULL, ...) {
  if (length(rank) > 1) {
    opar <- graphics::par(mfrow = grDevices::n2mfrow(length(rank)))
    graphics::par(mar = c(2, 2, 2, 1))
    on.exit(graphics::par(opar))
  }
  idx <- rank
  if (any(idx > length(x))) { stop("object shorter than requested numer of ranks") }
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
    # browser()
    # the legend could be established based on he summary result...
    # leg <- summary(x)[idx[i],]
    # leg_lab <- colnames(leg)
    # leg_val <- leg
    leg_lab <- colnames(attr(x[[idx[i]]], "scores"))
    #leg_val <- as.character(round(x = as.numeric(attr(x[[idx[i]]], "scores")), digits = c(3, 0, 3, 0, 3, 2, 2, 2, 2, 2)))
    leg_val <- sapply(1:10, function(j) { formatC(x = as.numeric(attr(x[[idx[i]]], "scores"))[j], digits = c(4, 0, 4, 0, 2, 2, 2, 2, 2, 2)[j], format = "f") })
    leg_val[2] <- rownames(attr(x[[idx[i]]], "scores"))
    graphics::legend("topright", legend = paste(leg_lab, leg_val, sep=": "), bty = "n", cex = 0.75, text.col = legend_text_col[, i])
    mhyp <- attr(x[[idx[i]]], "scores")[, "neutral_mass"]
    color <- if (!is.null(correct_mass)) {
      if (abs(mhyp - correct_mass) < 0.01) 3 else 2
    } else {
      1
    }
    graphics::mtext(sprintf("[%d] %.4f", idx[i], round(mhyp, 4)), col = color, cex = 0.9, line = 0.05)
  }
}

#'@rdname findMAIN
#'@param x Object of class findMAIN.
#'@param ... Further parameters.
#'@export
#'@method print findMAIN
print.findMAIN <- function (x, ...) {
  nres <- length(x)
  i <- 1
  scores <- summary(x)
  scores_i <- attr(x[[i]], "scores")
  nprec <- length(unique(sapply(x, function(x) attr(x, "scores")[,"adductmz"])))
  nadduct <- length(unique(sapply(x, function(x) attr(x, "scores")[,"adducthyp"])))
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
      "Selected m/z=%.4f as %s adduct of neutral mass %.4f with score %.2f.",
      scores_i[, 1],
      rownames(scores_i)[1],
      scores_i[, 3],
      scores_i[, 10]
    )
  )
  print(x[[i]])
}