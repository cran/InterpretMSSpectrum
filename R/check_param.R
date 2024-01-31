#' @title check_param.
#' @description Utility functions to check parameters.
#' @param x Character of length=1 providing a correct_peak info.
#' @param isotopes dataframe of allowed isotopes.
#' @param silent Suppress potential warning message.
#' @return A correct version of the parameter or NULL.
#' @keywords internal
#' @noRd
#' @importFrom enviPat check_chemform mergeform
#' @examples
#' inp <- "Name, C6H12O6, 0"
#' InterpretMSSpectrum:::check_param_correct_peak(x = inp)
#' # return NULL in case that input formula contains other elements than allowed by isotpes parameter
#' InterpretMSSpectrum:::check_param_correct_peak(x = inp, matrix("C", ncol=1))
#' InterpretMSSpectrum:::check_param_correct_peak(x = inp, matrix("C", ncol=1), silent=FALSE)
#' InterpretMSSpectrum:::check_param_correct_peak(x = "wrong, string", silent=FALSE)
#' InterpretMSSpectrum:::check_param_correct_peak(x = "wrong, string, really", silent=FALSE)
#'
check_correct_peak <- function(x = NULL, isotopes = NULL, silent = TRUE) {
  if (!is.null(x)) {
    test <- FALSE
    msg <- NULL
    if (is.null(isotopes)) isotopes <- as.matrix(InterpretMSSpectrum::chemical_elements[,1])
    tmp <- strsplit(x, ", ")[[1]]
    if (length(tmp) < 3) {
      test <- TRUE
      msg <- "can not split string into 3 components"
    }
    if (length(tmp) > 3) {
      tmp <- c(paste(tmp[1:(length(tmp) - 2)], collapse = ","), tmp[-(1:(length(tmp) - 2))])
    }
    if (!test && suppressWarnings(is.na(as.numeric(tmp[3])))) {
      test <- TRUE
      msg <- ifelse(is.null(msg), "can not convert 3rd component to numeric", paste(msg, "can not convert 3rd component to numeric", sep=", "))
    }
    fml_chk <- enviPat::check_chemform(isotopes = isotopes, chemforms = tmp[2])
    if (fml_chk$warning) {
      test <- TRUE
      msg <- ifelse(is.null(msg), "formula check of 2nd component failed", paste(msg, "formula check of 2nd component failed", sep=", "))
    } else {
      tmp[2] <- fml_chk$new_formula
      x <- paste(tmp, collapse = ", ")
      # store formula for later use separately
      attr(x, "mph") <- enviPat::mergeform(tmp[2], "H1")
    }
    if (test) {
      x <- NULL
      if (!silent) warning("InterpretMSSpectrum: correct_peak set to NULL as it did not pass QC (", msg, ").")
    }
  }
  return(x)
}

#' @name update_local_check
#' @keywords internal
#' @noRd
update_local_check <- function(local_check, correct_peak, rdisop_res, nval=0) {
  if (
    !is.null(correct_peak) && 
    local_check==0 && 
    length(grep(attr(correct_peak, "mph"), rdisop_res[[length(rdisop_res)]][,1])) != 1
  ) {
    local_check <- nval
  } 
  return(local_check)
}

#' @name check_param
#' @keywords internal
#' @noRd
check_param <- function(x = NULL) {
  
  param <- x
  
  # param.default <- list("ionization"="ESI",
  #                       "ionmode"="positive",
  #                       "allowed_elements"=c("C","H","N","O","P","S","Cl","Na","K"), # for Rdisop formula generation
  #                       "maxElements"="P4S4Na2K2", # for Rdisop formula generation
  #                       "minElements"="C1",
  #                       "substitutions"=data.frame("s1"=c("H1","H1","Na1","Na1","K1"),"s2"=c("Na1","K1","K1","H1","H1")),
  #                       "quick_isos"=TRUE,
  #                       "score_cutoff"=0.5,
  #                       "neutral_loss_cutoff"=2,
  #                       "ruleset"="none")
  # save(param.default, file = "Data/param.default.rda", compress = "gzip")
  
  stopifnot(is.list(param) | param %in% c("APCIpos", "ESIpos", "ESIneg", "default"))
  
  # load default parameters
  # using a parameter set is an attempt to summarize a number of parameters useful for either GC-APCI or LC-ESI
  # however, to be more flexible this parameter set can be provided directly as a list
  param.default <- InterpretMSSpectrum::param.default

  # now modify values of the default parameter set according to the provided option in 'param'
  if (is.list(param)) {
    for (n in names(param)) param.default[[n]] <- param[[n]]
  } else {
    if (param == "APCIpos") {
      param.default$"ionization" <- "APCI"
      param.default$"allowed_elements" <- c("C", "H", "N", "O", "P", "S", "Si")
      param.default$"maxElements" <- "P2S2"
      param.default["substitutions"] <- list(NULL)
      param.default$"neutral_loss_cutoff" <- 0.5
      param.default$"quick_isos" <- FALSE
      param.default$"ruleset" <- "APCI"
    }
    if (substr(param, 1, 3) == "ESI") {
      param.default$"allowed_elements" <- c("C", "H", "N", "O", "P", "S", "Cl")
      param.default$"maxElements" <- "P4S4"
      param.default$"substitutions" <- data.frame("s1" = c("H1", "H1", "Na1"), "s2" = c("Na1", "K1", "K1"))
      param.default$"ruleset" <- "ESI"
    }
    if (substr(param, nchar(param) - 2, nchar(param)) == "pos") {
      param.default$"ionmode" <- "positive"
    } else {
      param.default$"ionmode" <- "negative"
    }
  }
  param <- param.default
  param$"isotopes" <- as.matrix(param$"allowed_elements", ncol = 1)
  param$"iso_mass" <- ifelse(param$"ionization" == "APCI", 1.0015, 1.003355)
  
  return(param)
  
}

#' @name check_met_db
#' @keywords internal
#' @noRd
check_met_db <- function(x = NULL, isotopes = NULL, silent = TRUE) {
  met_db <- x
  if (is.null(isotopes)) isotopes <- as.matrix(InterpretMSSpectrum::chemical_elements[,1])
  if (!is.null(met_db)) {
    test <- FALSE
    colN <- grep("[Nn]ame", colnames(met_db))[1]
    if (length(colN) == 1) colnames(met_db)[colN] <- "Name"
    colF <- grep("[Ff]ormula", colnames(met_db))[1]
    if (length(colF) == 1) colnames(met_db)[colF] <- "Formula"
    colM <- which(colnames(met_db) %in% c("M+H", "mz"))[1]
    if (length(colM) == 1) colnames(met_db)[colM] <- "M+H"
    if (!all(c("Name", "Formula", "M+H") %in% colnames(met_db))) {
      test <- TRUE
    } else {
      met_db <- met_db[, c("Name", "Formula", "M+H"), drop = FALSE]
      met_db[, "Name"] <- as.character(met_db[, "Name"])
      met_db[, "Formula"] <- as.character(met_db[, "Formula"])
      met_db[, "M+H"] <- as.numeric(met_db[, "M+H"])
      met_db <- met_db[!is.na(met_db[, "M+H"]), , drop = FALSE]
      if (nrow(met_db) == 0) {
        test <- TRUE 
      } else {
        # remove rows with empty Formula 
        met_db <- met_db[!(is.na(met_db[, "Formula"]) | met_db[, "Formula"] == ""), , drop = FALSE]
      }
      if (nrow(met_db) == 0) {
        test <- TRUE
      } else {
        # remove rows with invalid Formula 
        met_db <- met_db[!enviPat::check_chemform(isotopes = isotopes, met_db[, "Formula"])$warning, , drop = FALSE]
      }
      if (nrow(met_db) == 0) test <- TRUE
    }
    if (test) {
      met_db <- NULL
      if (!silent) warning("InterpretMSSpectrum: met_db set to NULL as it didn't pass QC.")
    }
  }
  return(met_db)
}

#' @name check_sf_db
#' @examples
#' sf_db <- system.file("extdata", "APCI_min.db", package = "InterpretMSSpectrum")
#' check_sf_db(x = NULL)
#' check_sf_db(x = sf_db)
#' @keywords internal
#' @noRd
check_sf_db <- function(x = NULL, isotopes = NULL, silent = TRUE) {
  formula_db <- x
  if (!is.null(formula_db)) {
    if (is.character(formula_db) && length(formula_db) == 1 && file.exists(formula_db)) {
      check_pkg <- sapply(c("RSQLite", "DBI"), requireNamespace, quietly = TRUE)
      if (!all(check_pkg)) {
        msg <- paste0(
          "The use of a predefined Database requires package", ifelse(sum(!check_pkg) > 1, "s", ""),
          paste(names(check_pkg)[!check_pkg], collapse = ", "),
          ". Please install. Running with 'formula_db'=NULL."
        )
        formula_db <- NULL
      } else {
        # $ToDo$ We need to extend the isotopes list in the param object to include all
        # elements of formulas in the DB or alternatively remove some entries from DB before proceeding
        db_con <- DBI::dbConnect(RSQLite::SQLite(), formula_db)
        if (inherits(db_con, "SQLiteConnection")) {
          formula_db <- db_con
        } else {
          formula_db <- NULL
        }
      }
    }
  }
  return(formula_db)
}
  
#' @name check_neutral_losses
#' @examples
#' check_neutral_losses()
#' check_neutral_losses(NULL, "ESI")
#' @keywords internal
#' @noRd
check_neutral_losses <- function(x = NULL, ionization = c("APCI", "ESI"), isotopes = NULL, silent = TRUE) {
  ionization <- match.arg(ionization)
  typical_losses_definition <- x
  if (is.null(isotopes)) isotopes <- as.matrix(InterpretMSSpectrum::chemical_elements[,1])
  # load neutral loss table
  if (is.null(typical_losses_definition)) {
    if (ionization == "APCI") {
      neutral_losses <- InterpretMSSpectrum::neutral_losses_APCI
    }
    if (ionization == "ESI") {
      neutral_losses <- InterpretMSSpectrum::neutral_losses_ESI
    }
  } else {
    if (length(typical_losses_definition) == 1 && is.character(typical_losses_definition) && file.exists(typical_losses_definition)) {
      neutral_losses <- utils::read.table(typical_losses_definition, sep = "\t", header = T, dec = ",", as.is = T)
    } else {
      neutral_losses <- typical_losses_definition
    }
  }
  # ensure proper formula in neutral_loss table (for later add/sub-molecule functions)
  neutral_losses[, "Formula"] <- enviPat::check_chemform(isotopes = isotopes, neutral_losses[, "Formula"])[, "new_formula"]
  return(neutral_losses)
}
