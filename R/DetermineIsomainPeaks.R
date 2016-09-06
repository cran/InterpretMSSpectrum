#'@title DetermineIsomainPeaks.
#'
#'@description
#'\code{DetermineIsomainPeaks} will evaluate a mass spectrum and try to find the main isotopic clusters.
#'
#'@details
#'This function is used within \link{PlotSpec} and \link{InterpretMSSpectrum} to identify main isotopic clusters.
#'It is currently exported to allow the user to modify/substitute the result but may become an internal function in the future.
#'
#'@param spec A two-column matrix with mz and int.
#'@param int_cutoff Do not consider isomain peaks below this intensity cutoff (relative to base peak).
#'@param dmz_cutoff Expected maximum within scan mass defect of your device in milli Dalton.
#'@param precursor Specify the assumed precursor explicitly (ensure that precursor mass is included in list and everything above/higher is removed).
#'
#'@return
#'A vector of ion masses from a spectrum which are potential fragment masses (without isotopes).
#'
#'@keywords internal
#'
DetermineIsomainPeaks <-
function(spec=NULL, int_cutoff=0.03, dmz_cutoff=0.001, precursor=NULL) {

  # extract isotope groups by mass gap search
  isomain <- split(data.frame(spec), GetGroupFactor(round(spec[,1]), gap=1.1))

  # test for +H2O - CH4 shift (equals a +1.979265 shift which for some metabolites is more intense than the M+H)
  # and take otherwise maximum intensity mass per group (assuming that M+H is favored/of higher intensity relativ to M+)
  isomain <- sapply(isomain, function(x) { 
    max_int <- which.max(x[,2])
    test <- abs(x[max_int,1]-1.979265-x[,1]) <= dmz_cutoff
    test <- test & (x[,2]/max(x[,2]))>0.5
    ifelse(any(test), x[which(test),1], x[max_int,1])
  })

  # apply intensity cutoff for small peaks (may loose informative peaks but will speed up process)
  isomain <- which(spec[,1] %in% isomain & spec[,2] > int_cutoff*max(spec[,2]))
  
  if (is.null(precursor)) {
    # if no precursor is set, make some additional quality checks for isomain peaks
    
    # remove high mz peaks if their intensity <15% of base peak; empirical threshold for 'contaminating' high masses
    # !! modified in v05
    i_max <- isomain[which.max(spec[isomain,2])]
    while(spec[isomain[length(isomain)],2] < 0.15*spec[i_max,2]) {
      isomain <- isomain[-length(isomain)]
    }
    
    # remove TMS adducts (highest isomain has a lower equivalent with dmz=72.03953)
    test <- abs(spec[isomain[length(isomain)],1]-72.03953-spec[isomain,1]) <= dmz_cutoff
    if (any(test)) isomain <- isomain[-length(isomain)]
    
    # remove +H2O - H2 peak which can be present and large compared to M+H
    test <- abs(spec[isomain[length(isomain)],1]-15.99491-spec[isomain,1]) <= dmz_cutoff
    if (any(test)) isomain <- isomain[-length(isomain)]
    
    # test for (i) high mass (ii) high Int (iii) natural isotopes (iv) neutral loss partner
    # and limit number of total isomain-peaks (speed, memory)
    # ...
    #n <- 5
    #if (length(isomain)>n) { isomain <- isomain[order(spec[isomain,2],decreasing = TRUE)[1:n]] }
    
  } else {
    # ensure that precursor mass is included in list and everything above/higher is removed
    test <- which(abs(spec[isomain,1]-precursor) <= 1)
    if (length(test)==0) {
      if (any(abs(spec[,1]-precursor) <= 1)) {
        isomain <- c(isomain, which.min(abs(spec[,1]-precursor)))
        test <- length(isomain)
      }
    }
    if (length(test)>=2) {
      test <- test[which.min(abs(spec[isomain[test],1]-precursor))]
    }
    if (length(test)==1) {
      isomain <- isomain[!(spec[isomain,1] > spec[isomain[test],1])]
    }
  }  
  return(spec[isomain,1])
}