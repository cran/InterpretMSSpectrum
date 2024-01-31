#' @title Interpreting High-Res-MS spectra.
#'
#' @description \code{InterpretMSSpectrum} will read, evaluate and plot a deconvoluted
#'     mass spectrum (mass*intensity pairs) from either TMS-derivatized GC-APCI-MS data
#'     or ESI+/- data. The main purpose is to identify the causal metabolite or more
#'     precisely the sum formula of the molecular peak by annotating and interpreting
#'     all visible fragments and isotopes.
#'
#' @details For further details refer to and if using please cite
#'     Jaeger et al. (2016, <doi:10.1021/acs.analchem.6b02743>) in case of GC-APCI and
#'     Jaeger et al. (2017, <doi:10.1002/rcm.7905>) for ESI data. The Interpretation is
#'     extremely speed up if 'formula_db' (a predetermined database of potential sum
#'     formulas) is provided within the function call. Within the package you may use
#'     \link{GenerateMetaboliteSQLiteDB} to prepare one for yourself or request
#'     a download link from \email{jan.lisec@bam.de} as de-novo calculation for a wide
#'     mass range may take several days.
#'
#' @param spec A 2-column matrix of mz/int pairs. If spec=NULL then \code{InterpretMSSpectrum}
#'     tries to read data from clipboard (i.e. two columns copied from an Excel spreadsheet).
#' @param precursor The ion (m/z) from spec closest to this mass will be considered as
#'     precursor (can be nominal, i.e. if precursor=364 then 364.1234 would be selected from
#'     spectrum if it is closest).
#' @param correct_peak For testing purposes. A character in the form of "name, formula, mz"
#'     to evaluate spectra against. Note! Separating character is ', '.
#' @param met_db A metabolite DB (e.g. GMD or internal) can be provided to search for
#'     candidates comparing M+H ions (cf. Examples).
#' @param typical_losses_definition A file name (e.g. D:/BuildingBlocks_GCAPCI.txt) from
#'     where to load relevant neutral losses (cf. Details). Alternatively an data frame with
#'     columns 'Name', 'Formula' and 'Mass'.
#' @param silent Logical. If TRUE no plot is generated and no output except final candidate
#'     list is returned.
#' @param dppm Specifies ppm error for Rdisop formula calculation.
#' @param param Either a list of parameters or a character shortcut indicating a predefined
#'     set. Currently 'APCIpos', 'ESIpos' and 'ESIneg' are supported as shortcuts. If a list
#'     is provided, list elements will overwrite similar named entries from the list that
#'     can be accessed using utils::data(param.default, package = "InterpretMSSpectrum").
#' @param formula_db A pre calculated database of sum formulas and their isotopic fine
#'     structures can be used to extremely speed up the function.
#'
#' @return An annotated plot of the mass spectrum and detailed information within the
#'     console. Main result, list of final candidate formulas and their putative fragments,
#'     will be returned invisibly.
#'
#' @examples
#' # load APCI test data
#' apci_spectrum <- InterpretMSSpectrum::apci_spectrum
#'
#' # (otional) provide information of a correct peak as a character containing
#' # name, formula and ion mass -- separated by ',' as shown below
#' cp <- "Glutamic acid (3TMS), C14H33NO4Si3, 364.1790"
#'
#' # (otional) provide a database of known peaks
#' mdb <- data.frame(
#'   "Name" = c("Glutamic acid (3TMS)", "other peak with same sum formula"),
#'   "Formula" = c("C14H33NO4Si3", "C14H33NO4Si3"),
#'   "M+H" = c(364.179, 364.179), stringsAsFactors = FALSE, check.names = FALSE
#' )
#'
#' # (otional) provide a database of precalculated formulas to speed up the process
#' fdb <- system.file("extdata", "APCI_min.db", package = "InterpretMSSpectrum")
#'
#' # apply function providing above arguments which will print to the console
#' # and open a new plot
#' InterpretMSSpectrum(spec = apci_spectrum, correct_peak = cp, met_db = mdb, formula_db = fdb)
#'
#' @export
#'
#' @importFrom enviPat check_chemform mergeform

InterpretMSSpectrum <- function(
    spec = NULL,
    precursor = NULL,
    correct_peak = NULL,
    met_db = NULL,
    typical_losses_definition = NULL,
    silent = FALSE,
    dppm = 3,
    param = "APCIpos",
    formula_db = NULL) {
  # check for Rdisop to be able to keep it in suggested packages
  if (!requireNamespace("Rdisop", quietly = TRUE)) {
    message("\nThe use of this function requires package 'Rdisop'.\nPlease install with 'BiocManager::install(\"Rdisop\")\'\n")
    return(NULL)
  }

  # POTENTIAL PARAMETERS that could be allowed for the user to modify
  # !!! carefull limiting (to maximum of 5 peaks) is experimental
  max_isomain_peaks <- 5

  # Input CHECKS ----
  # check provided parameter list
  param <- check_param(x = param)

  # correct_peak
  correct_peak <- check_correct_peak(x = correct_peak, isotopes = param$isotopes, silent = silent)
  
  # met_db
  met_db <- check_met_db(x = met_db, isotopes = param$isotopes, silent = silent)
  
  # data_base
  sf_db <- check_sf_db(x = formula_db, isotopes = NULL, silent = TRUE)
  # register the disconnect of SQLite connection upon function exit
  if (!is.null(sf_db)) { on.exit(DBI::dbDisconnect(sf_db)) }
  
  neutral_losses <- check_neutral_losses(x = typical_losses_definition, isotopes = param$isotopes, silent = silent)
  
  # Internal FUNCTIONS ----
  GetFragmentData <- function(M0 = NULL, spec = NULL, n = 2, iso_mass = 1.003355) {
    # try to get reasonable isotope peaks for M0 from spectrum
    p <- sapply(0:n, function(dmz) {
      tmp <- which(abs((M0 + dmz * iso_mass) - spec[, 1]) < 0.01)
      if (length(tmp) > 1) {
        tmp <- which.min(abs((M0 + dmz * iso_mass) - spec[, 1]))
      }
      return(ifelse(length(tmp) == 1, tmp, NA))
    })
    if (any(is.na(p))) p <- p[1:(min(which(is.na(p))) - 1)]
    frag <- rbind(spec[p, 1], spec[p, 2] / sum(spec[p, 2], na.rm = T))
    attr(frag, "M0") <- M0
    return(frag)
  }
  EvaluateFragment <- function(frag = NULL, dppm = 3, param = NULL, silent = FALSE, sf_db = NULL) {
    # we need to correct our observed fragment data with the mass of an electron (charge) being +/- 0.00055 depending on positive/negative mode
    em <- 0.00055
    frag[1, ] <- frag[1, ] + ifelse(param$ionmode == "positive", 1, -1) * em
    # if (round(frag[1,1])==364) browser()
    if (!is.null(sf_db)) {
      # data base approach --> much faster compared to default solution below
      dmz <- 0.0005 + frag[1, 1] * dppm / 10^6
      dbq <- DBI::dbSendQuery(sf_db, "SELECT * FROM sfdb WHERE Mass > (:x) AND Mass < (:y);", data.frame(x = frag[1, 1] - dmz, y = frag[1, 1] + dmz))
      out <- DBI::dbFetch(dbq, -1)
      DBI::dbClearResult(dbq)
      if (nrow(out) >= 1) {
        # adaptive algorithm to compute isotopic pattern on the fly and fill sf_db over time...
        miss_iso <- is.na(out[, "m0"])
        if (any(miss_iso)) {
          out[miss_iso, 4:11] <- plyr::ldply(out[miss_iso, "Formula"], function(fml) {
            matrix(t(GetIsotopeDistribution(fml = fml, res = 30000, n = 3, ele_vec = param$allowed_elements)), nrow = 1)
          })
          # sf_db[mcand[miss_iso],4:11] <<- out[miss_iso,4:11]
        }
        # we predetermined 3 isotopes so we can compare at max 3 even if we find more within our spectrum...
        max_iso <- min(c(ncol(frag), 3))
        out[, "Score"] <- sapply(1:nrow(out), function(j) {
          the <- matrix(unlist(c(out[j, 4:(4 + max_iso - 1)], out[j, 8:(8 + max_iso - 1)])), nrow = 2, byrow = TRUE)
          the[2, ] <- the[2, ] / sum(the[2, ])
          mScore(obs = frag[, 1:max_iso, drop = FALSE], the = the, dppm = dppm)
        })
        out <- out[order(out[, "Score"], decreasing = TRUE), c("Formula", "Score", "Valid", "Mass"), drop = F]
      } else {
        out <- data.frame("Formula" = I(character(0)), "Score" = numeric(0), "Valid" = character(0), "Mass" = numeric(0))
      }
      attr(out, "M0") <- attr(frag, "M0")
    } else {
      # de novo calculation
      if (length(frag[1, ]) >= 2) {
        molecules <- Rdisop::decomposeIsotopes(frag[1, ], frag[2, ], mzabs = 0.0005, ppm = dppm, z = 1, maxisotopes = 1 + length(frag[1, ]), elements = Rdisop::initializeElements(param$allowed_elements), minElements = param$minElements, maxElements = param$maxElements)
      } else {
        molecules <- Rdisop::decomposeMass(frag[1, ], mzabs = 0.0005, ppm = dppm, z = 1, maxisotopes = 1, elements = Rdisop::initializeElements(param$allowed_elements), minElements = param$minElements, maxElements = param$maxElements)
      }
      if (is.null(molecules)) {
        out <- data.frame("Formula" = I(character(0)), "Score" = numeric(0), "Valid" = character(0), "Mass" = numeric(0))
        attr(out, "M0") <- attr(frag, "M0")
      } else {
        out <- data.frame("Formula" = I(molecules$formula), "Score" = molecules$score, "Valid" = molecules$valid, "Mass" = molecules$exactmass)
        attr(out, "M0") <- attr(frag, "M0")
        # extract isotope fine structure
        if (param$ionization == "APCI" & !param$quick_isos) {
          # for APCI data the Rdisop isotope distribution is less correct and can be better estimated using enviPat
          isos <- lapply(1:nrow(out), function(j) {
            x <- GetIsotopeDistribution(fml = out[j, "Formula"], res = 30000, n = 2, ele_vec = param$allowed_elements)
            # fall back sollution for no C suggestions
            if (ncol(x) <= 2) {
              x <- molecules[["isotopes"]][[j]]
            }
            if (diff(range(x[1, ])) <= (ncol(x) - 1.5)) {
              x <- molecules[["isotopes"]][[j]]
            }
            return(x)
          })
        }
        if (param$ionization == "ESI" | param$quick_isos) {
          isos <- lapply(molecules[["isotopes"]], function(x) {
            x[1, ] <- x[1, ]
            x[2, ] <- x[2, ] / sum(x[2, ])
            return(x)
          })
        }
        # the whole following part tries to determined an optimal and fair integrative score for each suggestion, including mass precision and intensity distribution for up to 2 isotopes
        # but it is slow and for some fraction of the spectra error prone...
        out[, "Score"] <- sapply(1:nrow(out), function(j) {
          max_iso <- min(c(ncol(frag), ncol(isos[[j]])))
          mScore(obs = frag[, 1:max_iso, drop = FALSE], the = isos[[j]][, 1:max_iso, drop = FALSE], dppm = dppm)
        })
        # ensure with enviPat valid chemical formulas (necessary for neutral loss detection and scoring later)
        out[, "Formula"] <- enviPat::check_chemform(isotopes = param$isotopes, out[, "Formula"])[, "new_formula"]
      }
    }
    out <- out[order(out[, "Score"], decreasing = TRUE), , drop = F]
    return(out)
  }
  RemoveEmptyFragments <- function(rdisop_res, silent = TRUE, step = "") {
    # remove fragments without suggestions
    flt <- sapply(rdisop_res, nrow) == 0
    if (any(flt)) {
      if (all(flt)) {
        warning(paste0("[RemoveEmptyFragments] No Fragments left after step ", step))
        # keep empty list of length=1
        rdisop_res <- rdisop_res[1]
      } else {
        if (!silent) paste0("[RemoveEmptyFragments] Fragments", paste(which(flt), collapse = " & "), "removed at step ", step)
        flt <- rev(which(sapply(rdisop_res, nrow) == 0))
        for (k in flt) rdisop_res[[k]] <- NULL
      }
    }
    return(rdisop_res)
  }
  GetRdisopResult <- function(spec = NULL, isomain = NULL, silent = TRUE, dppm = 3, param = NULL, sf_db = sf_db) {
    if (!is.null(attr(isomain, "adducthyp"))) {
      rdisop_res <- lapply(isomain, function(M0) {
        frag <- GetFragmentData(M0 = M0, spec = spec, n = 3, iso_mass = param$iso_mass)
        if (M0 == isomain[length(isomain)]) {
          # to calculate de nove sum formulas we need to correct the MID for the (de)protonation
          frag[1, ] <- frag[1, ] - attr(isomain, "adducthyp") + ifelse(param$ionmode == "positive", 1, -1) * 1.0073
        }
        rdisop_res <- EvaluateFragment(frag = frag, dppm = 2 * dppm, silent = silent, param = param, sf_db = sf_db)
        invisible(rdisop_res)
      })
    } else {
      rdisop_res <- lapply(isomain, function(M0) {
        frag <- GetFragmentData(M0 = M0, spec = spec, n = 3, iso_mass = param$iso_mass)
        rdisop_res <- EvaluateFragment(frag = frag, dppm = 2 * dppm, silent = silent, param = param, sf_db = sf_db)
        invisible(rdisop_res)
      })
    }
    RemoveEmptyFragments(rdisop_res, silent = silent, step = "GetRdisopResult")
  }
  RemoveByScore <- function(rdisop_res, score_cutoff = 0, silent = TRUE) {
    rdisop_res <- lapply(rdisop_res, function(x) {
      x[x[, "Score"] >= (score_cutoff * max(x[, "Score"])), ]
    })
    rdisop_res <- RemoveEmptyFragments(rdisop_res, silent = silent, step = "RemoveByScore")
  }
  GenerateMainOutput <- function(rdisop_res_list, stats, met_db, isomain) {
    cat(paste("\n\nTotal number of formulas per fragment before and after filtering...\n"))
    print(stats)
    rdisop_res_best <- rdisop_res_list[[1]]
    rdisop_res_best <- rdisop_res_best[!is.na(rdisop_res_best[, 1]), ]
    cat(paste("\n\nDetails of best candidate...\n"))
    print(rdisop_res_best)
    MH <- ifelse(nrow(rdisop_res_best) == 1, 1, ifelse(abs(rdisop_res_best[nrow(rdisop_res_best), "Mass"] - rdisop_res_best[nrow(rdisop_res_best) - 1, "Mass"] - 72.0395) < 0.0005, nrow(rdisop_res_best) - 1, nrow(rdisop_res_best)))
    # check resulting candidate list against a DB for best matche
    if (!is.null(met_db)) {
      met_db[, "Formula"] <- enviPat::check_chemform(isotopes = param$isotopes, met_db[, "Formula"])[, "new_formula"]
      M0 <- enviPat::subform(rdisop_res_best[MH, "Formula"], "H1")
      if (any(met_db$Formula %in% M0)) {
        if (!silent) print(met_db[met_db$Formula %in% M0, ])
        best_cand <- paste(met_db[met_db$Formula %in% M0, "Name"], collapse = "; ")
        best_cand_col <- 3
      } else {
        if (any(abs(met_db$"M+H" - rdisop_res_best[MH, "Mass"]) < 0.002)) {
          if (!silent) print(met_db[abs(met_db$"M+H" - rdisop_res_best[MH, "Mass"]) < 0.002, -1])
          best_cand <- paste(met_db[abs(met_db$"M+H" - rdisop_res_best[MH, "Mass"]) < 0.002, "Name"], collapse = "; ")
          best_cand_col <- 6
        } else {
          best_cand <- ""
          best_cand_col <- grDevices::grey(0.9)
        }
      }
    } else {
      best_cand <- ""
      best_cand_col <- 1
    }
    # outcommented to allow combination of spectrum plot with other plots
    # opar <- par(no.readonly=TRUE)
    # set colors for spectra plot (isomain peaks get a red color)
    tmp.col <- rep(1, nrow(spec))
    tmp.col[which(spec[, 1] %in% isomain)] <- 2
    fml_expr <- sapply(rdisop_res_best[, "Formula"], function(x) {
      x <- CountChemicalElements(x)
      x[x==1] <- ""
      paste0(names(x),"[", x, "]", collapse="*")
    },USE.NAMES = FALSE)
    fml_expr <- gsub("[[]]", "", fml_expr)
    PlotSpec(x = spec, cols = tmp.col, txt = data.frame("mz" = rdisop_res_best[, "Mass"], "Formula" = fml_expr, "expr"=TRUE), neutral_losses = neutral_losses, neutral_loss_cutoff = param$neutral_loss_cutoff, substitutions = param$substitutions, ionization = param$ionization)
    graphics::mtext(paste("Remaining combinations:", length(rdisop_res_list)), line = -1.2, adj = 0, side = 3, col = grDevices::grey(0.5))
    graphics::mtext(best_cand, line = -2.4, adj = 0, side = 3, col = best_cand_col)
    if (!is.null(correct_peak)) {
      graphics::mtext(correct_peak, line = -3.6, adj = 0, side = 3, col = grDevices::grey(0.5))
      fcor <- strsplit(correct_peak, ", ")[[1]]
      fcor <- fcor[-1][grep("^C[[:digit:]]", fcor[-1])]
      if (length(fcor) == 1) {
        fcor <- enviPat::check_chemform(isotopes = param$isotopes, chemforms = fcor)
        if (fcor[, 1]) {
          warning("[InterpretMSSpectrum] Probably a wrong specification of 'correct_peak'", call. = FALSE)
        } else {
          fcor <- enviPat::mergeform(fcor[, 2], "H1")
        }
        # [Modification:if correct peak is specified but smaller than the 'observed' M0 the now uncommented line leads to trouble]
        # fcor <- which(sapply(rdisop_res_list, function(x) {fcor %in% x[,"Formula"]}))
        fcor <- which(sapply(rdisop_res_list, function(x) {
          fcor == x[nrow(x), "Formula"]
        }))
        cat(paste("\n\nRank of specified M0 =", ifelse(length(fcor) == 1, fcor, "NA"), "\n"))
        graphics::mtext(paste0("Rank = ", fcor), line = -4.6, adj = 0, side = 3, col = grDevices::grey(0.5))
      }
    }
    # restore parameters
    # par(opar)
  }

  # Main PART ----
  # read data (mz,int -table) from clipboard if not provided explicitly and sort by mz (should be standard but you never know...)
  if (is.null(spec)) spec <- ReadSpecClipboard()
  spec <- spec[order(spec[, 1]), , drop = FALSE]

  # check  if spectrum is okay
  if (!nrow(spec) >= 1) {
    # [ToDo] Further checks may be useful
    warning("[InterpretMSSpectrum] Spectra does not contain any information (nrow=0).", call. = FALSE)
    invisible(NULL)
  } else {
    # keep global error message for testing purposes (if correct_peak is known/provided)
    if (!is.null(correct_peak)) {
      # global_err <<- NULL # global assignment removed to pass RCheck without NOTEs
      global_err <- NULL
      local_check <- 0
    }

    # evaluate main peaks
    isomain <- DetermineIsomainPeaks(spec = spec, int_cutoff = 0.03, precursor = precursor, ionization = param$ionization, ionmode = param$ionmode, limit = max_isomain_peaks)

    # # modify spectra if no good [M+H] is found but alternative adduct hypothesis instead
    # if (!is.null(attr(isomain,"adducthyp"))) {
    #
    # }

    stats <- data.frame("mz" = round(isomain, 4), "initial" = NA, "score_cutoff" = NA, "PlausibleFormula" = NA, "TypicalLosses" = NA)

    # start timing for testing purposes
    time_elapse <- Sys.time()
    # use Rdisop to get potential formulas (using up to n=2 isotopic peaks found)
    rdisop_res <- GetRdisopResult(spec = spec, isomain = isomain, silent = silent, dppm = dppm, param = param, sf_db = sf_db)
    stats[stats[, "mz"] %in% round(sapply(rdisop_res, attr, "M0"), 4), "initial"] <- sapply(rdisop_res, nrow)
    local_check <- update_local_check(local_check = 0, correct_peak, rdisop_res, nval = 1)
    time_elapse <- c(time_elapse, Sys.time())

    # remove according to individual score based on mz deviation and isotopic fit
    rdisop_res <- RemoveByScore(rdisop_res, score_cutoff = param$score_cutoff, silent = silent)
    stats[stats[, "mz"] %in% round(sapply(rdisop_res, attr, "M0"), 4), 3] <- sapply(rdisop_res, nrow)

    local_check <- update_local_check(local_check, correct_peak, rdisop_res, nval = 2)
    time_elapse <- c(time_elapse, Sys.time())

    # restrict to plausible formulas
    if (is.null(sf_db)) {
      # if a predetermined database is used there is no need to check formula plausibility again
      rdisop_res <- lapply(rdisop_res, function(x) {
        x[sapply(x[, "Formula"], PlausibleFormula, ruleset = param$ruleset), ]
      })
      rdisop_res <- RemoveEmptyFragments(rdisop_res, silent = silent, step = "PlausibleFormula")
      stats[stats[, "mz"] %in% round(sapply(rdisop_res, attr, "M0"), 4), 4] <- sapply(rdisop_res, nrow)
    } else {
      # stats[stats[,"mz"] %in% round(sapply(rdisop_res,attr,"M0"),4),4] <- sapply(rdisop_res,nrow)
    }

    local_check <- update_local_check(local_check, correct_peak, rdisop_res, nval = 3)
    time_elapse <- c(time_elapse, Sys.time())

    # restrict based on neutral losses (only if more than 1 fragment is left in spectrum)
    # test of all losses are potentially helpful
    nl_vec <- sapply(neutral_losses[, "Formula"], function(x) {
      neutral_losses[neutral_losses[, "Formula"] == x, "Mass"]
    })
    if (length(rdisop_res) >= 2) {
      for (k in 1:length(nl_vec)) {
        rdisop_res <- RestrictByTypicalLosses(rdisop_res = rdisop_res, tl = nl_vec[k], neutral_loss_cutoff = param$neutral_loss_cutoff, punish = 0.5, substitutions = param$substitutions)
        local_check <- update_local_check(local_check, correct_peak, rdisop_res, nval = 4)
      }
    }
    stats[stats[, "mz"] %in% round(sapply(rdisop_res, attr, "M0"), 4), 5] <- sapply(rdisop_res, nrow)
    time_elapse <- c(time_elapse, Sys.time())

    if (sum(sapply(rdisop_res, nrow)) > 0) {
      # obtain most likely combination
      rdisop_res_list <- ScoreFormulaCombination(rdisop_res, nl_vec = nl_vec, punish_invalid = 0.5, punish_S = 0.2, punish_nonplausible = 0.5, return_rank = NA, neutral_loss_cutoff = param$neutral_loss_cutoff, substitutions = param$substitutions, silent = silent)
      time_elapse <- c(time_elapse, Sys.time())

      # plot annotated spectrum
      if (!silent) GenerateMainOutput(rdisop_res_list, stats, met_db, isomain)
      time_elapse <- c(time_elapse, Sys.time())

      # write potential error source to global
      if (!is.null(correct_peak) && local_check > 0) {
        global_err <- paste(correct_peak, "not present after", c("initial generation", "score filter", "plausibility filter", "neutral loss filter")[local_check])
        print(global_err)
      }

      time_elapse <- round(diff(time_elapse), 4)
      names(time_elapse) <- c("FormulaGen", "ScoreFilt", "Plausible", "NeutralLoss", "PathEval", "Plot")
      if (!silent) {
        cat("\n\nTime elapsed during individual processing steps...\n")
        print(time_elapse)
      }

      # return final result list
      attr(rdisop_res_list, "stats") <- stats
      invisible(rdisop_res_list)
    } else {
      invisible(NULL)
    }
  }
}
