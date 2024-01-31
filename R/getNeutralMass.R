#' @title getNeutralMass.
#' @description \code{getNeutralMass} will return the neutral mass calculated for a measured m/z and a specified potential adduct.
#' @details To achieve the task formulas are split into elements and counted using \link{CountChemicalElements}.
#' @param m Measured m/z.
#' @param adduct_name Name of potential adduct.
#' @param adduct_rules data frame of possible adducts.
#' @return Numeric vector.
#'
#' @examples
#' InterpretMSSpectrum:::getNeutralMass()
#' rules <- data.frame(
#'   "name" = c("[M+1]+","[2M+1]+","[M+1]2+","[2M+1]2+","[M+2]+","[2M+2]+","[M+2]2+","[2M+2]2+"), 
#'   "nmol" = c(1,2,1,2,1,2,1,2), 
#'   "charge" = c(1,1,2,2,1,1,2,2), 
#'   "massdiff" = c(1,1,1,1,2,2,2,2)
#' )
#' nm <- sapply(1:nrow(rules), function(i) { 
#'   InterpretMSSpectrum:::getNeutralMass(100, rules[i,1], rules)
#' })
#' cbind(mz=100, nm, rules)
#' 
#' # set up testing rules for negaive ionization mode
#' rules_neg <- rules
#' rules_neg$name <- gsub("[+]", "-", rules_neg$name)
#' rules_neg$charge <- -1 * rules_neg$charge
#' rules_neg$massdiff <- -1 * rules_neg$massdiff
#' nm <- sapply(1:nrow(rules_neg), function(i) { 
#'   InterpretMSSpectrum:::getNeutralMass(100, rules_neg[i,1], rules_neg)
#' })
#' cbind(mz=100, nm, rules_neg)
#' 
#'
#' @keywords internal
#' @noRd
#'
getNeutralMass <- function(
  m=100, 
  adduct_name = "[2M-2H]2-", 
  adduct_rules=data.frame("name"="[2M-2H]2-", "nmol"=2, charge=2, massdiff=-2.014553)) 
{
  i <- which(adduct_rules[,"name"]==adduct_name)
  if (any(i)) {
    nmol <- adduct_rules[i,"nmol"]
    charge <- abs(adduct_rules[i,"charge"])
    massdiff <- adduct_rules[i,"massdiff"]
  } else {
    nmol <- 1
    charge <- 1
    massdiff <- 0
  }
  return(charge*m/nmol-massdiff)
}