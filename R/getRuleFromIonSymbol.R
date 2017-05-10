#' Generate adduct rule from ion symbol
#' 
#' Translate an ion symbol to an adduct rule. This function is used internally by \code{findMAIN}, but may be useful elsewhere.
#' 
#' @param ions character vector of ion symbols, e.g. "[M+H]+", "[M+Na]+", 
#'   "[M+H-NH3]-". Please use full notation in square brackets, though some 
#'   frequent ions can be abbrevated ("M+H","M+Na","M+K","M+NH4", "M+", "M", 
#'   "M-H","M+Cl-", "M-").
#'   
#' @return A data frame with four columns "name", "nmol", "charge", "massdiff".
#' @export
#' @keywords internal
#'   
#' @examples
#' getRuleFromIonSymbol(c("[M+H]+", "[M+Na]+"))
getRuleFromIonSymbol <- function(ions="[M+H]+") {
    checkSymbol <- function(ion) {
        regexpr("\\[[0-9]{0,2}M.*\\][0-9]{0,2}[\\+\\-]{1,2}", ion) != -1
    }
    shortCuts <- cbind(c("M+H","M+Na","M+K","M+NH4", "M+", "M",
                         "M-H","M+Cl-", "M-"),
                       c("[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+", "[M]+", "[M]+",
                         "[M-H]-","[M+Cl]-", "[M]-"))
    em <- 0.0005485799
    chemical_elements <- NULL
    utils::data(chemical_elements, envir=environment(), package="InterpretMSSpectrum")
    on.exit(rm(chemical_elements))
    out <- lapply(ions, function(ion) {
      if(ion %in% shortCuts[,1]) ion <- shortCuts[,2][ which(shortCuts[,1] == ion) ]
      if(!checkSymbol(ion)) stop("invalid ion")
      nmol <- sub(".*[^0-9M]([0-9]?M).*", "\\1", ion)
      nmol <- sub("M", "", nmol)
      nmol <- as.numeric(ifelse(nmol=="", 1, nmol))
      ch <- sub(".*[^0-9]([0-9]{0,2}[\\+\\-])$", "\\1", ion)
      sgn <- sub("[^\\+\\-]", "", ch)
      sgn <- ifelse(sgn=="+", 1, -1)
      ch <- sub("[\\+\\-]", "", ch)
      ch[ch==""] <- "1"
      ch <- as.numeric(ch)
      ch <- ch * sgn
      x <- ion
      x <- sub("^.*\\[", "", x)
      x <- sub("\\].*", "", x)
      x <- sub("[0-9]?M", "", x)
      starts <- gregexpr("[\\+\\-]", x)[[1]]
      ends <- c(starts[-1]-1, nchar(x))
      n <- length(starts)
      spl <- lapply(1:n, function(i) substr(x, starts[i], ends[i]))
      massdiff <- lapply(spl, function(y) {
        sgn <- sub("^([\\+\\-]).*", "\\1", y)
        sgn <- ifelse(sgn=="+", 1, -1)
        el <- sub("^[\\+\\-]", "", y)
        if (regexpr("^[0-9]+[A-Za-z]+", el) != -1) el <- gsub("([0-9]+)([A-Za-z]+)", "\\2\\1", el)
        #browser()
        #el <- tabulateElements(el)
        el <- CountChemicalElements(x=el)
        masses <- sapply(names(el), function(a) {
          chemical_elements[,2][ which(chemical_elements[,1] == a)[1] ]
        })
        return(sum(masses * el) * sgn)
      })
      massdiff <- sum(unlist(massdiff), na.rm=TRUE) + ch * -em
      return(data.frame(name=ion, nmol=nmol, charge=ch, massdiff=massdiff, stringsAsFactors = FALSE))
    })
    return(do.call("rbind", out))
}
