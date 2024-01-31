#' @title get_exactmass.
#' 
#' @description Get the exact mass for chemical formulas.
#' 
#' @param x Vector of chemical formulas.
#'   
#' @return A named vector of exact masses.
#' 
#' @export
#'   
#' @examples
#' InterpretMSSpectrum::get_exactmass(c("C6H12O6", "Na", "H1"))
#' 
get_exactmass <- function(x) {
  ce <- InterpretMSSpectrum::chemical_elements
  ce <- stats::setNames(ce$mass,ce$name)
  sapply(x, function(fml) {
    fml <- CountChemicalElements(x = fml)
    return(sum(ce[names(fml)]*fml))
  })
}