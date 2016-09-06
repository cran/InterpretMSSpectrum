#'@title GetGroupFactor.
#'
#'@description
#'\code{GetGroupFactor} will provide a factor vector which allows to split any numeric vector
#' into Groups according to a gap parameter for further processing.
#'
#'@details
#'not exported
#'
#'@param x numeric vector.
#'@param gap difference up from which a new group is assumed.
#'
#'@return
#'factor of length x, with groups of x which are seperated by a value >gap.
#'
#'@keywords internal
#'
GetGroupFactor <-
function(x, gap) {
  stopifnot(is.numeric(x))
  idx <- rank(x)
  x <- x[order(x)]
	x <- c(T, diff(x)>gap)
  x <- factor(rep(1:sum(x), times=diff(c(which(x),length(x)+1))))
	return(x[idx])
}