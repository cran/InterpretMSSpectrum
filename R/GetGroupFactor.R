#'@title GetGroupFactor.
#'
#'@description
#'\code{GetGroupFactor} will split a numeric vector according to a specified gap value. This is often a useful tool and therefore exported to the namespace.
#'
#'@param x Numeric vector.
#'@param gap Difference between two consecutive values at which a split is generated.
#'
#'@return
#'A factor vector of length(x) indicating the diffferent groups in x.
#'
#'@examples
#'x <- c(1:3,14:12,6:9)
#'GetGroupFactor(x=x, gap=2)
#'split(x, GetGroupFactor(x=x, gap=2))
#'
#'@export
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