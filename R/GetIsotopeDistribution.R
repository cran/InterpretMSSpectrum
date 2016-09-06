#'@title GetIsotopeDistribution.
#'
#'@description
#'\code{GetIsotopeDistribution} will generate an isotopic distribution for a given formula.
#'
#'@details
#'not exported
#'
#'@param fml sum formula.
#'@param res MS resolution. Yet experimental, may fail.
#'@param n Number of isotopes to calculate.
#'@param ele_vec Character vector of elements to consider.
#'
#'@return
#'Isotope distribution formatted similar to Rdisop result but more precise using enviPat.
#'
#'@import enviPat
#'@importFrom stats weighted.mean
#'@importFrom utils data
#'
#'@keywords internal
#'
GetIsotopeDistribution <- function(fml=NULL, res=NULL, n=2, ele_vec=c("C","H","N","O","P","S","Si")) {
  # load and restrict isotope list locally
  utils::data("isotopes", package="enviPat", envir=environment())
  isotopes <- isotopes[as.character(isotopes[,"element"]) %in% ele_vec & isotopes[,"abundance"]>=0.001,]
  # check formula
  fml <- enviPat::check_chemform(isotopes, chemforms=fml)$new_formula
  # calculate and transform isotopic pattern
  if (is.null(res)) {
    isopat <- enviPat::isopattern(isotopes = isotopes, chemforms = fml, threshold=0, verbose = FALSE)[[1]]
    g <- GetGroupFactor(x=isopat[,1], gap=0.2)
    theo <- sapply(levels(g), function(x) { c(round(stats::weighted.mean(x = isopat[g==x,1], w = isopat[g==x,2]),4), sum(isopat[g==x,2]/100)) })
  } else {
    isopat <- enviPat::isopattern(isotopes = isotopes, chemforms = fml, threshold=0.0001, verbose=FALSE)
    env <- enviPat::envelope(isopat, resolution=30000, verbose = FALSE)
    ipt <- enviPat::vdetect(env, detect="intensoid", plotit=FALSE, verbose = FALSE)
    theo <- t(ipt[[1]])
  }
  theo <- theo[,1:min(c(ncol(theo),(n+1))),drop=F]
  theo[2,] <- round(theo[2,]/sum(theo[2,]),4)
  return(theo)
}