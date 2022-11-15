#' @title is.subformula.
#'
#' @description \code{is.subformula} will test for all elements of one chemical formula (f_sub) to be present in another (f_main).
#'
#' @details To achieve the task formulas are split into elements and counted using \link{CountChemicalElements}.
#'
#' @param f_sub Supposed chemical sub formula.
#' @param f_main Supposed chemical main formula.
#' @param substitutions data frame of allowed substitutions to consider.
#'
#' @return Logical indicating if 'f_sub' is potentially a sub formula of 'f_main'.
#'
#' @examples
#' InterpretMSSpectrum:::is.subformula(f_sub = "C6H12O6", f_main = "C6H12O6")
#' InterpretMSSpectrum:::is.subformula(f_sub = "C4H8O5", f_main = "C6H12O6")
#'
#' @keywords internal
#' @noRd
#'
is.subformula <- function(f_sub, f_main, substitutions=NULL) {
  ExtrChemElem <-  function (x) {
    p <- gregexpr("[[:upper:]]", x)[[1]]
    return(unique(gsub("[^[:alpha:]]", "", sapply(1:length(p), function(i) { substr(x, p[i], ifelse(i == length(p), nchar(x), p[i + 1] - 1)) }))))
  }
  stopifnot(is.character(f_sub), is.character(f_main))
  ele <- unique(c(ExtrChemElem(f_sub), ExtrChemElem(f_main)))
  s <- CountChemicalElements(x=f_sub, ele = ele)
  m <- CountChemicalElements(x=f_main, ele = ele)
  if (!is.null(substitutions)) {
    count_subs <- lapply(1:nrow(substitutions), function(i) { lapply(substitutions[i,],CountChemicalElements) })
    count_subs <- count_subs[sapply(count_subs, function(x) { all(names(x[[1]])%in%names(s)) && s[names(x[[1]])]>=x[[1]] && all(names(x[[2]])%in%names(m)) && m[names(x[[2]])]>=x[[2]] })]
    if (length(count_subs)>=1) {
      test <- sapply(count_subs, function(x) {
        s[names(x[[1]])] <- s[names(x[[1]])]-x[[1]]
        m[names(x[[2]])] <- m[names(x[[2]])]-x[[2]]
        all(m >= s)
      })
      if (any(test)) {
        x <- count_subs[[which(test)[1]]]
        s[names(x[[1]])] <- s[names(x[[1]])]-x[[1]]
        m[names(x[[2]])] <- m[names(x[[2]])]-x[[2]]
      }
    }
  }
  out <- all(m >= s)
  names(out) <- paste(names(m)[(m-s)>0], (m-s)[(m-s)>0], sep="", collapse="")
  return(out)
}