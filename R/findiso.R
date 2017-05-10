#'@title findiso.
#'
#'@description
#'\code{findiso} will evaluate a mass spectrum and try to find the main isotopic clusters and adducts.
#'
#'@details
#'This function is used within \link{findMAIN} to identify main isotopic clusters.
#'
#'@param spec A two-column matrix with mz and int.
#'@param intthr Do not consider isomain peaks below this intensity cutoff (relative to base peak).
#'@param mzabs Expected maximum within scan mass defect of your device in Dalton.
#'@param CAMERAlike T/F reformatting the result similar to a CAMERA isocluster search output.
#'
#'@return
#'An annotated spectrum.
#'
#'@keywords internal
#'
findiso <- function(spec=NULL, mzabs=0.001, intthr=0.03, CAMERAlike=TRUE) {
  ## helper function for gap search (usually integrated in InterpretMS)
  GetGroupFactor <-
    function(x, gap) {
      stopifnot(is.numeric(x))
      idx <- rank(x)
      x <- x[order(x)]
      x <- c(T, diff(x)>gap)
      x <- factor(rep(1:sum(x), times=diff(c(which(x),length(x)+1))))
      return(x[idx])
    }
  ## helper function to determine main mz and charge of an isotopic group
  ## [may still need work]
  mainmz <- function(x, l=1) {
    ## l (=loading) not yet used
    max_int <- which.max(x[,2])
    test <- abs(x[max_int,1]-1.979265-x[,1]) <= mzabs
    test <- test & (x[,2]/max(x[,2]))>0.5
    if(any(test)) {
      x[which(test),1]
    } else {
      if (max_int>=4) {
        x[which.max(x[1:4,2]),1]
      } else {
        x[max_int,1]
      }
    }
  }
  s <- as.data.frame(spec)
  s[,"type"] <- NA
  ## extract isotope groups by mass gap search
  isomain <- split(s, GetGroupFactor(round(s[,1]), gap=1.1))
  ## limit relative to base peak
  isomain <- isomain[sapply(isomain,function(x){any(x[,2]>(intthr*max(s[,2])))})]
  ## check multiple loading and annotate main peaks
  for (i in 1:length(isomain)) {
    x <- isomain[[i]]
    x <- x[ x[,2] > max(x[,2])/1000, , drop=FALSE]
    if (nrow(x)>=2) {
      for (l in rev(1:3)) {
        filt <- sapply(x[,1], function(y) { any(abs(x[,1]-y-1.003355/l)<0.01 | abs(x[,1]-y+1.003355/l)<0.01) })
        if (any(filt)) {
          m <- x[filt,,drop=FALSE]
          ## keep only those with relation to biggest peak
          mmax <- which.max(m[,2])
          m <- m[sapply(m[,1], function(y) { any(abs(m[mmax,1]-y-seq(-5,5)*1.003355/l)<0.01) }),]
          m0 <- mainmz(x=m)
          ## restrict to determined main and its isotopes
          n <- max(round(m[,"mz"]-m0))
          m <- m[sapply(m0+(0:n)*(1.003355/l), function(mz) { which.min(abs(m[,"mz"]-mz))}),]
          s[s[,1]==m0,"type"] <- paste0(c("M","D","T")[l],nrow(m)-1)
          x <- x[!filt | 1:nrow(x)%in%max(which(filt)),,drop=FALSE]
        }
      }
    } else {
      s[s[,1]==x[,1],"type"] <- "S"
    }
  }
  ## reformat to match CAMERA annotation scheme
  if (CAMERAlike) {
    s[,"isogr"] <- NA
    s[,"iso"] <- NA
    s[,"charge"] <- NA
    mainiso <- which(!(is.na(s[,"type"]) | s[,"type"]=="S"))
    if (any(mainiso)) {
      for (i in 1:length(mainiso)) {
        ty <- s[mainiso[i],"type"]
        ch <- (1:3)[c("M","D","T") %in% substr(ty,1,1)]
        n <- as.numeric(substr(ty,2,3))
        if (n>=1) {
          iso <- c(mainiso[i], sapply(s[mainiso[i],"mz"]+(1:n)*(1.003355/ch), function(mz){ which.min(abs(s[,"mz"]-mz))}))
        } else {
          iso <- mainiso[i]
        }
        s[iso,"isogr"] <- i
        s[iso,"iso"] <- (1:length(iso))-1
        s[iso,"charge"] <- ch
      }
    }
    s <- s[,!(colnames(s)%in%"type")]
  }
  ## keep correct peak attribute
  if (!is.null(attr(spec, "correct_peak"))) {
    attr(s, "correct_peak") <- attr(spec, "correct_peak")
    ##print(attr(s, "correct_peak"))
  }
  return(s)
}