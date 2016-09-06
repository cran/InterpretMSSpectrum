#'@title is.subformula.
#'
#'@description
#'\code{is.subformula} will test for all elements of one chemical formula (f_sub) to be present in another (f_main).
#'
#'@details
#'To achieve the task formulas are split into elements and counted using \link{CountChemicalElements}.
#'
#'@param f_sub Supposed chemical sub formula.
#'@param f_main Supposed chemical main formula.
#'
#'@return
#'Logical indicating if f_sub is potentially a subformula of f_main.
#'
#'@keywords internal
#'
is.subformula <- function(f_sub, f_main) {
  stopifnot(is.character(f_sub), is.character(f_main))
  count_f_sub <- CountChemicalElements(x=f_sub)
  count_f_main <- CountChemicalElements(x=f_main, ele=names(count_f_sub))
  return(all(count_f_main >= count_f_sub))
}