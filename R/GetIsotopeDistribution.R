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
#'@param check.fml The fml needs to be in enviPat style, i.e. not CH4 but C1H4. This will be ensured but can be skipped setting check to FALSE to speed up.
#'@param vdetect.detect Will be deprecated soon. Only for testing this enviPat parameter.
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
GetIsotopeDistribution <- function(fml=NULL, res=NULL, n=2, ele_vec=c("C","H","N","O","P","S","Si"), check.fml=TRUE, vdetect.detect=c("centroid","intensoid")[1]) {
  # load and restrict isotope list locally
  utils::data("isotopes", package="enviPat", envir=environment())
  isotopes <- isotopes[as.character(isotopes[,"element"]) %in% ele_vec & isotopes[,"abundance"]>=0.001,]
  # ensure formula to be in enviPat style
  fml <- enviPat::check_chemform(isotopes, chemforms=fml)$new_formula
  # calculate and transform isotopic pattern
  if (is.null(res)) {
    isopat <- enviPat::isopattern(isotopes = isotopes, chemforms = fml, threshold=0, verbose = FALSE)[[1]]
    g <- GetGroupFactor(x=isopat[,1], gap=0.2)
    theo <- sapply(levels(g), function(x) { c(round(stats::weighted.mean(x = isopat[g==x,1], w = isopat[g==x,2]),4), sum(isopat[g==x,2]/100)) })
  } else {
    isopat <- enviPat::isopattern(isotopes = isotopes, chemforms = fml, threshold=0, verbose=FALSE)
    env <- enviPat::envelope(isopat, resolution=res, verbose = FALSE)
    ipt <- enviPat::vdetect(env, detect=vdetect.detect, plotit=FALSE, verbose = FALSE)
    theo <- t(ipt[[1]])
    #browser()
    theo <- sapply(theo[1,1]+(0:n)*1.003, function(mz){ theo[,which.max(abs(theo[1,]-mz)<0.1)] })
  }
  theo <- theo[,1:min(c(ncol(theo),(n+1))),drop=F]
  theo[2,] <- round(theo[2,]/sum(theo[2,]),4)
  return(theo)
}