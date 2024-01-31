#' @title InterpretTP.
#'
#' @description \code{InterpretTP} is a wrapper function around 
#'     \link{InterpretMSSpectrum} which will read, evaluate and plot 
#'     a deconvoluted mass spectrum (mass*intensity pairs) from either 
#'     TMS-derivatized GC-APCI-MS data or ESI+/- data. 
#'     It allows to provide a chemical formula as a potential precursor of the 
#'     spectrum. This formula will be used to set the parameters 
#'     'allowed_elements' and 'maxElements' during de-novo formula generation.
#'
#' @details For further details refer to \link{InterpretMSSpectrum}.
#'
#' @inheritDotParams InterpretMSSpectrum
#' @param fml A chemical formula of the standard used for transformation product generation.
#' @param param Keyword or parameter list, similar as in \link{InterpretMSSpectrum}.
#'
#' @return An annotated plot of the mass spectrum and detailed information 
#'     within the console. Main result, list of final candidate formulas and 
#'     their putative fragments, will be returned invisibly.
#'
#' @examples
#' # load test data
#' utils::data(apci_spectrum)
#'
#' # provide information of a correct peak (if you know) as character
#' cp <- "Glutamic acid (3TMS), C14H33NO4Si3, 364.1790"
#'
#' # provide database of known peaks and correct peak
#' mdb <- data.frame(
#'   "Name" = c("Glutamic acid (3TMS)", "other peak with same sum formula"),
#'   "Formula" = c("C14H33NO4Si3", "C14H33NO4Si3"),
#'   "M+H" = c(364.179, 364.179), stringsAsFactors = FALSE, check.names = FALSE
#' )
#' 
#' # provide a database of precalculated formulas to speed up the process
#' fdb <- system.file("extdata", "APCI_min.db", package = "InterpretMSSpectrum")
#'
#' # apply function providing above arguments (dppm is set to 0.5 to reduce run time)
#' InterpretTP(fml = "C14H33NO4Si3", spec=apci_spectrum, param="APCIpos")
#'
#' @export
#'
#' @importFrom enviPat check_chemform mergeform

InterpretTP <- function(fml=NULL, param="APCIpos", ...) {
  
  stopifnot(is.vector(fml) && is.character(fml) && length(fml)==1)
  stopifnot(is.list(param) | param %in% c("APCIpos","ESIpos","ESIneg","default"))

  # load default parameters
  # using a parameter set is an attempt to summarize a number of parameters useful for either GC-APCI or LC-ESI
  # however, to be more flexible this parameter set can be provided directly as a list
  param.default <- InterpretMSSpectrum::param.default
  
  ele <- CountChemicalElements(x = fml)
  
  # now modify values of the default parameter set according to the provided option in 'param'
  if (is.list(param)) {
    for (n in names(param)) param.default[[n]] <- param[[n]]
  } else {
    if (substr(param,nchar(param)-2,nchar(param))=="pos") {
      param.default$"ionmode" <- "positive"
    } else {
      param.default$"ionmode" <- "negative"
    }
    if (param=="APCIpos") {
      param.default$"ionization" <- "APCI"
      if ("H" %in% names(ele)) ele["H"] <- ele["H"]+1 else ele <- c(ele, "H"=1)
      param.default["substitutions"] <- list(NULL)
      param.default$"neutral_loss_cutoff" <- 0.5
      param.default$"quick_isos" <- FALSE
      param.default$"ruleset" <- "APCI"
    }
    if (substr(param,1,3)=="ESI") {
      if (param.default$"ionmode"=="positive") for (e in c("H", "Na", "K")) if (e %in% names(ele)) ele[e] <- ele[e]+1 else { ele <- c(ele, 1); names(ele)[length(ele)] <- e }
      param.default$"substitutions" <- data.frame("s1"=c("H1","H1","Na1"), "s2"=c("Na1","K1","K1"))
      param.default$"ruleset" <- "ESI"
    }
  }
  param.default$"allowed_elements" <- names(ele)
  param.default$"maxElements" <- paste0(names(ele), ele, collapse="")
  param <- param.default
  param$"isotopes" <- as.matrix(param$"allowed_elements", ncol=1)
  param$"iso_mass" <- ifelse(param$"ionization"=="APCI", 1.0015, 1.003355)

  return(InterpretMSSpectrum(param=param, ...))
  
}
