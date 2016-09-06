#'@title mScore.
#'
#'@description
#'\code{mScore} will calculate a mass defect weighted score for an mz/int values measure for an isotopic cluster in comparison to the theoretically expected pattern.
#'
#'@details
#'The maximum expected average mass error should be specified in ppm. A observed pattern deviating
#'that much from the theoretical pattern would still receive a reasonable (average) mScore while
#'observations deviating stronger or less strong will reach lower or higher mScores respectively.
#'Likewise the intensity presision should specify the average quality of your device to maintain
#'stable isotopic ratios.
#'
#'@param obs Observed (measured) values, a matrix with two rows (mz/int).
#'@param the Theoretical (estimated from sum formula) values, a matrix with two rows (mz/int).
#'@param mass_prec The expected mass precision will influence mScore (see Details).
#'@param int_prec The expected intensity precision will influence mScore (see Details).
#'@param limit minimal value of mScore. Should be left on zero.
#'@param rnd_prec Rounding precision of mScore.
#'
#'@return
#'Scalar mScore giving the quality of the observed data if theoretical data are true.
#'
#'@examples
#'# get theoretical isotopic pattern of Glucose
#'glc <- Rdisop::getMolecule("C6H12O6")$isotopes[[1]][,1:3]
#'mScore(obs=glc, the=glc)

#'# simulate mass and int defects
#'ef <- function(x, e) {runif(1,x-x*e,x+x*e)}
#'glc_obs <- glc
#'glc_obs[1,] <- sapply(glc[1,], ef, e=2*10^-6)
#'glc_obs[2,] <- sapply(glc[2,], ef, e=0.02)
#'mScore(obs=glc_obs, the=glc)

#'# simulate mass and int defects systematically
#'ef <- function(x, e) {runif(1,x-x*e,x+x*e)}
#'n <- 11
#'mz_err <- round(seq(0,5,length.out=n),3)
#'int_err <- round(seq(0,0.1,length.out=n),3)
#'mat <- matrix(NA, ncol=n, nrow=n, dimnames=list(mz_err, 100*int_err))
#'for (i in 1:n) {
#'  glc_obs[1,] <- sapply(glc[1,], ef, e=mz_err[i]*10^-6)
#'  for (j in 1:n) {
#'    glc_obs[2,] <- sapply(glc[2,], ef, e=int_err[j])
#'    mat[i,j] <- mScore(obs=glc_obs, the=glc)
#'  }
#'}
#'plot(x=1:n, y=1:n, type="n",axes=FALSE, xlab="mass error [ppm]", ylab="isoratio error [%]")
#'axis(3,at=1:n,rownames(mat),las=2); axis(4,at=1:n,colnames(mat),las=2); box()
#'cols <- grDevices::colorRampPalette(colors=c(2,6,3))(diff(range(mat))+1)
#'cols <- cols[mat-min(mat)+1]
#'text(x=rep(1:n,each=n), y=rep(1:n,times=n), labels=as.vector(mat), col=cols)
#'
#'@export
#'
mScore <- function(obs=NULL, the=NULL, mass_prec=2, int_prec=0.05, limit=0, rnd_prec=0) { 
  # the average mass/int quality of the data is specified within the function call
  # this is multiplied with a factor internally to render average data with mScore above 50
  qfac <- 2
  stopifnot(all(dim(obs)==dim(the)))
  max_err_mz <- qfac*mass_prec*the[1,]/10^6
  dmz <- 1+99*abs(obs[1,]-the[1,])/max_err_mz
  max_err_int <- qfac*int_prec*the[2,]
  dint <- 1+99*abs(obs[2,]-the[2,])/max_err_int
  out <- round(101-mean(sqrt(dmz*dint)), rnd_prec)
  return(ifelse(out<limit, limit, out))
}
