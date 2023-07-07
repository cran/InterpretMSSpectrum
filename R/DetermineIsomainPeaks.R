#' @title DetermineIsomainPeaks.
#'
#' @description
#' \code{DetermineIsomainPeaks} will evaluate a mass spectrum and try to find the main isotopic clusters.
#'
#' @details
#' This function is used within \link{PlotSpec} and \link{InterpretMSSpectrum} to identify main isotopic clusters.
#' It is currently exported to allow the user to modify/substitute the result but may become an internal function in the future.
#'
#' @param spec A two-column matrix with mz and int.
#' @param int_cutoff Do not consider isomain peaks below this intensity cutoff (relative to base peak).
#' @param dmz_cutoff Expected maximum within scan mass defect of your device in Dalton.
#' @param precursor Specify the assumed precursor explicitly (ensure that precursor mass is included in list and everything above/higher is removed).
#' @param ionization Should be made ionization dependent (not ready yet).
#' @param ionmode Either 'positive' or 'negative'.
#' @param limit Limit final list to a maximum number of peaks to speed up follow up processes.
#'
#' @return
#' A vector of ion masses from a spectrum which are potential fragment masses (without isotopes).
#'
#' @keywords internal
#' @noRd
#'
#' @examples
#' #load test data and apply function
#' \donttest{
#'   utils::data(apci_spectrum, package = "InterpretMSSpectrum")
#'   InterpretMSSpectrum:::DetermineIsomainPeaks(spec=apci_spectrum, ionization="APCI")
#'   utils::data(esi_spectrum, package = "InterpretMSSpectrum")
#'   InterpretMSSpectrum:::DetermineIsomainPeaks(spec=esi_spectrum, ionization="ESI")
#' }
DetermineIsomainPeaks <- function(spec=NULL, int_cutoff=0.03, dmz_cutoff=0.001, precursor=NULL, ionization=NULL, ionmode="positive", limit=NULL) {
  
  if (ionization=="ESI") {
    # analyze spectrum
    isomain <- findiso(spec=spec, mzabs=dmz_cutoff, intthr=int_cutoff, CAMERAlike=TRUE)
    # exclude multiple charged peaks and isotopes
    if (any(isomain[,"charge"]>=2,na.rm=T)) isomain <- isomain[-which(isomain[,"charge"]>=2),]
    if (any(isomain[,"iso"]>=1,na.rm=T)) isomain <- isomain[-which(isomain[,"iso"]>=1),]
    # if no precursor is set, try to guess it from data
    if (is.null(precursor)) {
      fmr <- summary(InterpretMSSpectrum::findMAIN(spec=spec, ionmode = ionmode))[1,,drop=FALSE]
      precursor <- fmr[,"adductmz"]
      adducthyp <- getRuleFromIonSymbol(fmr[,"adducthyp"])[1,4]
      attr(spec,"pot_mh") <- precursor-adducthyp+1.0073
      attr(spec,"adducthyp") <- adducthyp
      # exclude masses higher than mass which is closest to precursor (i.e. function is a bit fuzzy)
      MpH_pos <- which.min(abs(isomain[,"mz"]-attr(spec,"pot_mh")))
      prec_pos <- which.min(abs(isomain[,"mz"]-precursor))
      flt <- sort(c(which(isomain[,"mz"]<=isomain[MpH_pos,"mz"]),prec_pos))
    } else {
      flt <- which(isomain[,"mz"]<=isomain[which.min(abs(isomain[,"mz"]-precursor)),"mz"])
    }
    #browser()
    isomain <- isomain[flt,]
    n <- which(isomain[,"iso"]==0)
    # ensure that precursor is within list
    if (!nrow(isomain) %in% n) n <- c(n, nrow(isomain))
    # export this list according to limit settings
    if (is.null(limit)) limit <- length(isomain[,2])
    if (length(n)>limit) {
      # take at least all peaks with isotopes despite a limit parameter
      isomain <- isomain[n,"mz"]
    } else {
      #take high abundant peaks (assuming they have isotopes) and the precursor
      n <- order(isomain[,2], decreasing = TRUE)[1:limit]
      if (!nrow(isomain) %in% n) n <- c(n, nrow(isomain))
      isomain <- sort(isomain[n,"mz"])
    }
    if (length(isomain)>=1) {
      out <- isomain
    } else {
      #fallback sollution if nothing remained
      out <- spec[which.max(spec[,2]),1]
    }
  }
  
  if (ionization=="APCI") {
    # extract isotope groups by mass gap search
    isomain <- split(data.frame(spec), GetGroupFactor(round(spec[,1]), gap=1.1))
    
    # limit relative to base peak
    isomain <- isomain[sapply(isomain,function(x){any(x[,2]>(int_cutoff*max(spec[,2])))})]
    
    # exclude multiple charged peaks from groups
    # isomain <- lapply(isomain, function(x) {
    #   if (nrow(x)>=2) {
    #     filt <- sapply(x[,1], function(y) {!any(abs(x[,1]-y-0.5017)<0.01 | abs(x[,1]-y-0.3345)<0.01)})
    #     # remove also last isotope which could not be checked
    #     if (filt[length(filt)] & !filt[length(filt)-1]) filt[length(filt)] <- FALSE
    #     return(x[filt,,drop=FALSE])
    #   } else {
    #     return(x)
    #   }
    # })
    # isomain <- isomain[sapply(isomain,nrow)!=0]
    
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
    
      # remove high mz peaks if their intensity <15% of base peak; empirical threshold for 'contaminating' high masses
      # !! modified in v05
      i_max <- isomain[which.max(spec[isomain,2])]
      test <- spec[isomain[length(isomain)],2] < 0.15*spec[i_max,2]
      while(test) {
        # but don't do so if it has a CH4 loss [todo:] or they are connected by known neutral losses
        test2 <- any(abs(spec[isomain[length(isomain)],1]-16.0313-spec[isomain,1]) <= dmz_cutoff)
        if (test2) {
          test <- FALSE
        } else {
          isomain <- isomain[-length(isomain)]
          test <- spec[isomain[length(isomain)],2] < 0.15*spec[i_max,2]
        }
      }
      
      # remove TMS adducts (highest isomain has a lower equivalent with dmz=72.03953)
      test <- abs(spec[isomain[length(isomain)],1]-72.03953-spec[isomain,1]) <= dmz_cutoff
      if (any(test) && spec[isomain[test],2]>spec[isomain[length(isomain)],2]) isomain <- isomain[-length(isomain)]
      
      # remove +H2O - H2 peak which can be present and large compared to M+H
      test <- abs(spec[isomain[length(isomain)],1]-15.99491-spec[isomain,1]) <= dmz_cutoff
      if (any(test)) isomain <- isomain[-length(isomain)]
      
      # test for (i) high mass (ii) high Int (iii) natural isotopes (iv) neutral loss partner
  
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
    
    if (!is.null(limit) && length(isomain)>limit) {
      # limit number of total isomain-peaks (speed, memory)
      isomain <- sort(isomain[order(spec[isomain,2], decreasing=TRUE)[1:limit]])
      if (!any(spec[isomain,1]%in%precursor)) isomain <- c(isomain, which(spec[,1]%in%precursor))
    }
    out <- spec[isomain,1]
  }
  
  if (!is.null(attr(spec,"pot_mh"))) attr(out,"pot_mh") <- attr(spec,"pot_mh")
  if (!is.null(attr(spec,"adducthyp"))) attr(out,"adducthyp") <- attr(spec,"adducthyp")
  return(out)
  
}