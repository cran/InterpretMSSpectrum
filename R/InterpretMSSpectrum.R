#'@title Interpreting High-Res-MS spectra.
#'
#'@description
#'\code{InterpretMSSpectrum} will read, evaluate and plot a deconvoluted mass spectrum (mass*intensity pairs) from either TMS-derivatized GC-APCI-MS data or ESI+/- data. 
#'The main purpose is to identify the causal metabolite or more precisely the sum formula of the molecular peak by annotating and interpreting all visible fragments and isotopes.
#'
#'@details
#'For further details refer to and if using please cite Jaeger et al. (\url{http://dx.doi.org/10.1021/acs.analchem.6b02743}) in case of GC-APCI and Jaeger et al. 2017 (RCM, accepted) for ESI data.
#'The Interpretation is extremely speed up if 'formula_db' (a predetermined database of potential sum formulas) is provided within the function call. Within the package you may use \link{GenerateMetaboliteSQLiteDB} to prepare one for yourself or request a download link from \email{jan.lisec@charite.de} as de novo claculation may take several days.
#'
#'@param spec A 2-column matrix of mz/int pairs. If spec=NULL then \code{InterpretMSSpectrum} tries to read data from clipboard (i.e. two columns copied from an Excel speadsheet).
#'@param precursor The ion (m/z) from spec closest to this mass will be considered as precursor (can be nominal, i.e. if precursor=364 then 364.1234 would be selected from spectrum if it is closest).
#'@param correct_peak For testing purposes. A character in the form of "name, formula, mz" to evaluate spectra against. Note! Sperating character is ', '.
#'@param met_db A metabolite DB (e.g. GMD or internal) can be provided to search for candidates comparing M+H ions (cf. Examples).
#'@param typical_losses_definition A file name (e.g. D:/BuildingBlocks_GCAPCI.txt) from where to load relevant neutral losses (cf. Details). Alternatively an dataframe with columns Name, Formula and Mass.
#'@param silent Logical. If TRUE no plot is generated and no output except final candidate list is returned.
#'@param dppm Specifies ppm error for Rdisop formula calculation.
#'@param score_cutoff Specifies initial filtering step threshold per fragment. Sum Formulas with score_i < score_cutoff*max(score) will be removed.
#'@param neutral_loss_cutoff Specifies the allowed deviation in mDa for neutral losses to be accepted from the provided neutral loss list. If NULL determined dependent on ionization.
#'@param ionization Currently 'APCI' or 'ESI-' or 'ESI+' are supported as ionization modes.
#'@param quick_isos There are two ways to compute isotope pattern for comparison, a fast but error prone using Rdisop and a slow but (more) correct using enviPat.
#'@param formula_db A predetermined database of sum formulas and their isotopic fine structures can be used to extremely speed up the function.
#'
#'@return
#'An annotated plot of the mass spectrum and detailed information within the console.
#'Main result, list of final candidate formulas and their putative fragments, will be returned invisibly.
#'
#'@examples
#'#load test data
#'utils::data(apci_spectrum)
#'
#'# provide information of a correct peak (if you know)
#'correct_peak <- "Glutamic acid (3TMS), C14H33NO4Si3, 364.1790"
#'
#'# provide database of known peaks and correct peak
#'met_db <- data.frame("Name"=c("Glutamic acid (3TMS)","other peak with same sum formula"),
#'                     "Formula"=c("C14H33NO4Si3","C14H33NO4Si3"), 
#'                     "M+H"=c(364.179,364.179), stringsAsFactors=FALSE, check.names=FALSE)
#'                     
#'# apply function providing above arguments (dppm is set to 0.5 to reduce run time)
#'res <- InterpretMSSpectrum(spec=apci_spectrum, correct_peak=correct_peak, met_db=met_db, dppm=0.5)
#'
#'# show final function result (score-sorted list of potential fragment trees)
#'str(res)
#'\donttest{
#'#Given that you installed a prepared formula data base you can check performance increase by
#'setwd("D:/Bruker/R/Rpackage_InterpretMSSpectrum/")
#'system.time(InterpretMSSpectrum(spec=apci_spectrum, dppm=0.5, formula_db="APCI.db"))
#'system.time(InterpretMSSpectrum(spec=apci_spectrum, dppm=0.5, formula_db=NULL))
#'data(esi_spectrum)
#'plot(InterpretMSSpectrum::findMAIN(spec=esi_spectrum))
#'system.time(
#'InterpretMSSpectrum(spec=esi_spectrum, dppm=0.5, ionization="ESI+", formula_db="ESI.db")
#')
#'system.time(
#'InterpretMSSpectrum(spec=esi_spectrum, dppm=0.5, formula_db=NULL)
#')
#'}
#'
#'@export
#'
#'@import Rdisop
#'@import enviPat
#'@importFrom graphics mtext
#'@importFrom grDevices grey
#'@importFrom utils data read.table
#'
InterpretMSSpectrum <-
function(spec=NULL, precursor=NULL, correct_peak=NULL, met_db=NULL, typical_losses_definition=NULL, silent=FALSE, dppm=3, score_cutoff=0.5, neutral_loss_cutoff=NULL, ionization=c("APCI","ESI")[1], quick_isos=TRUE, formula_db=NULL) {

  # POTENTIAL PARAMETERS that could be allowed for the user to modify
  # !!! carefull limiting (to maximum of 5 peaks) is experimental
  max_isomain_peaks <- 5
  
  if (substr(ionization,nchar(ionization),nchar(ionization))%in%c("+","-")) {
    ionmode <- ifelse(substr(ionization,nchar(ionization),nchar(ionization))=="+","positive","negative")
    ionization <- substr(ionization,1,nchar(ionization)-1)
  } else {
    ionmode <- "positive"
  }

  if (ionization=="ESI") {
    #allowed_elements <- c("C","H","N","O","P","S","Cl","Na","K")
    #maxElements="P4S4Na2K2"
    #substitutions <- data.frame("s1"=c("H1","H1","Na1","Na1","K1"),"s2"=c("Na1","K1","K1","H1","H1"))
    allowed_elements <- c("C","H","N","O","P","S","Cl")
    maxElements="P4S4"
    substitutions <- data.frame("s1"=c("H1","H1","Na1"),"s2"=c("Na1","K1","K1"))
    if (is.null(neutral_loss_cutoff)) neutral_loss_cutoff <- 2
    #quick_isos <- TRUE
  }
  if (ionization=="APCI") {
    allowed_elements <- c("C","H","N","O","P","S","Si")
    maxElements="P2S2"
    substitutions <- NULL
    if (is.null(neutral_loss_cutoff)) neutral_loss_cutoff <- 0.5
    #quick_isos <- FALSE
  }
  em <- 0.00055
  #browser()
  
  # Input CHECKS
  # correct_peak
  if (!is.null(correct_peak)) {
    test <- FALSE
    #find_sep <- gregexpr(", ", correct_peak)[[1]]
    tmp <- strsplit(correct_peak,", ")[[1]]
    if (length(tmp)<3) test <- TRUE; msg <- "cant split string into 3 components"
    if (length(tmp)>3) tmp <- c(paste(tmp[1:(length(tmp)-2)],collapse=","),tmp[-(1:(length(tmp)-2))]) 
    if (!test && is.na(as.numeric(tmp[3]))) test <- TRUE; msg <- "cant convert 3rd component to numeric"
    #if (!test && enviPat::check_chemform(isotopes=data.frame("element"=allowed_elements,"isotope"="X","mass"=0,"abundance"=0), tmp[2])$warning) test <- TRUE
    if (!test && enviPat::check_chemform(isotopes=as.matrix(allowed_elements,ncol=1), tmp[2])$warning) {
      test <- TRUE
      msg <- "formula check of 2nd component failed"
    } else {
      tmp[2] <- enviPat::check_chemform(isotopes=as.matrix(allowed_elements,ncol=1), tmp[2])$new_formula
      # store formula for later use seperately
      fml <- enviPat::mergeform(tmp[2], "H1")
      correct_peak <- paste(tmp, collapse=", ")
    }
    if (test) {
      correct_peak <- NULL
      if (!silent) warning("InterpretMSSpectrum: correct_peak set to NULL as it didn't pass QC. (", msg, ")")
    }
  }
  # met_db
  if (!is.null(met_db)) {
    test <- FALSE
    colN <- grep("[Nn]ame", colnames(met_db))[1]; if (length(colN)==1) colnames(met_db)[colN] <- "Name"
    colF <- grep("[Ff]ormula", colnames(met_db))[1]; if (length(colF)==1) colnames(met_db)[colF] <- "Formula"
    colM <- which(colnames(met_db) %in% c("M+H","mz"))[1]; if (length(colM)==1) colnames(met_db)[colM] <- "M+H"
    if (!all(c("Name","Formula","M+H") %in% colnames(met_db))) {
      test <- TRUE
    } else {
      met_db <- met_db[,c("Name","Formula","M+H"),drop=FALSE]
      met_db[,"Name"] <- as.character(met_db[,"Name"])
      met_db[,"Formula"] <- as.character(met_db[,"Formula"])
      met_db[,"M+H"] <- as.numeric(met_db[,"M+H"])
      met_db <- met_db[!is.na(met_db[,"M+H"]),,drop=FALSE]
      if (nrow(met_db)==0) test <- TRUE else met_db <- met_db[!(is.na(met_db[,"Formula"]) | met_db[,"Formula"]==""),,drop=FALSE]
      if (nrow(met_db)==0) test <- TRUE else met_db <- met_db[!enviPat::check_chemform(isotopes=as.matrix(allowed_elements,ncol=1), met_db[,"Formula"])$warning,,drop=FALSE]
      if (nrow(met_db)==0) test <- TRUE
    }
    if (test) {
      met_db <- NULL
      if (!silent) warning("InterpretMSSpectrum: met_db set to NULL as it didn't pass QC.")
    }
  }
  # data_base
  sf_db <- FALSE
  if (!is.null(formula_db)) {
    if (is.character(formula_db) && length(formula_db)==1 && file.exists(formula_db)) {
      db_con <- DBI::dbConnect(RSQLite::SQLite(), formula_db)
      sf_db <- ifelse(class(db_con)=="SQLiteConnection", TRUE, FALSE)
    }
    # if (prod(dim(formula_db))>1 && all(c("Formula","Mass","m0","a0") %in% colnames(formula_db))) {
    #   # allow data.frame in parallel to SQLite DB ???
    # }
  }
  
  
  # load neutral loss table
  if (is.null(typical_losses_definition)) {
    #utils::data(paste("neutral_losses",ionization,sep="_"), envir=environment())
    if (ionization=="APCI") {
      neutral_losses_APCI <- NULL
      utils::data(neutral_losses_APCI, envir=environment(), package = "InterpretMSSpectrum")
      neutral_losses <- neutral_losses_APCI
    }
    if (ionization=="ESI") {
      neutral_losses_ESI <- NULL
      utils::data(neutral_losses_ESI, envir=environment(), package = "InterpretMSSpectrum")
      neutral_losses <- neutral_losses_ESI
    }
  } else {
    if (length(typical_losses_definition)==1 && is.character(typical_losses_definition) && file.exists(typical_losses_definition)) {
      neutral_losses <- utils::read.table(typical_losses_definition, sep="\t", header=T, dec=",", as.is=T)
    } else {
      neutral_losses <- typical_losses_definition
    }
  }
  # ensure proper formula in neutral_loss table (for later add/sub-molecule functions)
  neutral_losses[,"Formula"] <- enviPat::check_chemform(isotopes=as.matrix(allowed_elements,ncol=1), neutral_losses[,"Formula"])[,"new_formula"]
  
  # internal functions
  ReadSpecClipboard <- function() {
    # source could be Excel (German/English) or DA directly
    spec <- readLines("clipboard")
    spec <- gsub("\t"," ",spec) # replace Tabs
    if (length(grep("[^[:digit:],/. ]", spec[1]))==1) spec <- spec[-1] # strip header if present
    spec <- gsub(",",".",spec) # replace Colons
    spec <- gsub(" +$","",spec) # replace white space end
    spec <- gsub("^ +","",spec) # replace white space end
    spec <- t(as.matrix(sapply(spec, function(x) { as.numeric(strsplit(x," ")[[1]]) }))) # convert to numeric matrix
    if (ncol(spec)>=3) spec <- spec[,-1]
    return(spec)
  }
  GetFragmentData <- function(M0=NULL, spec=NULL, n=2, ionization="APCI") {
    # try to get reasonable isotope peaks for M0 from spectrum
    iso_mass <- ifelse(ionization=="APCI",1.001,1.003)
    p <- sapply(0:n, function(dmz) { 
      tmp <- which(abs((M0+dmz*iso_mass)-spec[,1])<0.01)
      if (length(tmp)>1) {
        tmp <- which.min(abs((M0+dmz*iso_mass)-spec[,1]))
      }
      return(ifelse(length(tmp)==1,tmp,NA)) 
    })
    if (any(is.na(p))) p <- p[1:(min(which(is.na(p)))-1)]
    frag <- rbind(spec[p,1], spec[p,2]/sum(spec[p,2], na.rm=T))
    attr(frag, "M0") <- M0
    return(frag)
  }
  EvaluateFragment <- function(frag=NULL, em=0.00055, dppm=3, score_cutoff=0, allowed_elements=c("C","H","N","O","P","S","Si"), ionization="APCI", silent=FALSE, maxElements="C60", quick_isos=FALSE, sf_db=NULL) {
    # we need to correct our observed fragment data with the mass of an electron (charge) being +/- 0.00055 depending on positive/negative mode
    frag[1,] <- frag[1,]+em
    #browser()
    if (sf_db) {
      # data base approach --> much faster compared to default sollution below
      dmz <- 0.0005+frag[1,1]*dppm/10^6
      dbq <- DBI::dbSendQuery(db_con, "SELECT * FROM sfdb WHERE Mass > (:x) AND Mass < (:y);", data.frame(x=frag[1,1]-dmz, y=frag[1,1]+dmz))
      out <- DBI::dbFetch(dbq, -1)
      DBI::dbClearResult(dbq)
      if (nrow(out)>=1) {
        # adaptive algorithm to compute isopattern on the fly and fill sf_db over time...
        miss_iso <- is.na(out[,"m0"])
        if (any(miss_iso)) {
          #browser()
          out[miss_iso,4:11] <- plyr::ldply(out[miss_iso,"Formula"], function(fml){matrix(t(GetIsotopeDistribution(fml=fml, res=30000, n=3, ele_vec=allowed_elements)),nrow=1)})
          #sf_db[mcand[miss_iso],4:11] <<- out[miss_iso,4:11]
        }
        # we predetermined 3 isotopes so we can compare at max 3 even if we find more within our spectrum...
        max_iso <- min(c(ncol(frag),3))
        out[,"Score"] <- sapply(1:nrow(out), function(j) {
          the <- matrix(unlist(c(out[j,4:(4+max_iso-1)],out[j,8:(8+max_iso-1)])),nrow=2,byrow=TRUE)
          the[2,] <- the[2,]/sum(the[2,])
          mScore(obs=frag[,1:max_iso,drop=FALSE], the=the, dppm=dppm)
        })
        out <- out[order(out[,"Score"], decreasing = TRUE),c("Formula","Score","Valid","Mass"),drop=F]
      } else {
        out <- data.frame("Formula"=I(character(0)),"Score"=numeric(0),"Valid"=character(0),"Mass"=numeric(0))
      }
      attr(out, "M0") <- attr(frag, "M0")
    } else {
      # de novo calculation
      if (length(frag[1,])>=2) {
        molecules <- Rdisop::decomposeIsotopes(frag[1,], frag[2,], mzabs=0.0005, ppm=dppm, z=1, maxisotopes=1+length(frag[1,]), elements=Rdisop::initializeElements(allowed_elements), minElements="C1", maxElements=maxElements)
      } else {
        molecules <- Rdisop::decomposeMass(frag[1,], mzabs=0.0005, ppm=dppm, z=1, maxisotopes=1, elements=Rdisop::initializeElements(allowed_elements), minElements="C1", maxElements=maxElements)
      }
      if (is.null(molecules)) {
        out <- data.frame("Formula"=I(character(0)),"Score"=numeric(0),"Valid"=character(0),"Mass"=numeric(0))
        attr(out, "M0") <- attr(frag, "M0")
      } else {
        out <- data.frame("Formula"=I(molecules$formula), "Score"=molecules$score, "Valid"=molecules$valid, "Mass"=molecules$exactmass)
        attr(out, "M0") <- attr(frag, "M0")
        # extract isotope fine structure
        if (ionization=="APCI" & !quick_isos) {
          # for APCI data the Rdisop isotope distribution is less correct and can be better estimated using enviPat
          isos <- lapply(1:nrow(out), function(j) { 
            x <- GetIsotopeDistribution(fml=out[j,"Formula"], res=30000, n=2, ele_vec=allowed_elements) 
            # fall back sollution for no C suggestions
            if (ncol(x)<=2) {
              x <- molecules[["isotopes"]][[j]]
            }
            if (diff(range(x[1,]))<=(ncol(x)-1.5)) {
              x <- molecules[["isotopes"]][[j]]
            }
            return(x)
          })
        }
        if (ionization=="ESI" | quick_isos) {
          isos <- lapply(molecules[["isotopes"]], function(x) {
            x[1,] <- x[1,]
            x[2,] <- x[2,]/sum(x[2,])
            return(x)
          })
        }
        # the whole following part tries to determined an optimal and fair integrative score for each suggestion, including mass precision and intensity distribution for up to 2 isotopes
        # but it is slow and for some fraction of the spectra error prone...
        out[,"Score"] <- sapply(1:nrow(out), function(j) {
          max_iso <- min(c(ncol(frag),ncol(isos[[j]])))
          mScore(obs=frag[,1:max_iso,drop=FALSE], the=isos[[j]][,1:max_iso,drop=FALSE], dppm=dppm)
        })
        # ensure with enviPat valid chemical formulas (necessary for neutral loss detection and scoring later)
        out[,"Formula"] <- enviPat::check_chemform(isotopes=as.matrix(allowed_elements,ncol=1), out[,"Formula"])[,"new_formula"]
      }
    }
    out <- out[order(out[,"Score"], decreasing = TRUE),,drop=F]
    return(out)
  }
  RemoveEmptyFragments <- function(rdisop_res, silent=TRUE, step="") {
    # remove fragments without suggestions
    flt <- sapply(rdisop_res,nrow)==0
    if (any(flt)) {
      if (all(flt)) {
        warning(paste0("[RemoveEmptyFragments] No Fragments left after step ", step))
        # keep empty list of length=1
        rdisop_res <- rdisop_res[1]
      } else {
        if (!silent) paste0("[RemoveEmptyFragments] Fragments", paste(which(flt),collapse=" & "), "removed at step ", step)
        flt <- rev(which(sapply(rdisop_res,nrow)==0))
        for (k in flt) rdisop_res[[k]] <- NULL
      }
    }
    return(rdisop_res)
  }
  GetRdisopResult <- function(spec=NULL, isomain=NULL, silent=TRUE, em=0.00055, dppm=3, allowed_elements=allowed_elements, ionization=ionization, maxElements=maxElements, quick_isos=quick_isos, sf_db=sf_db) {
    if (!is.null(attr(isomain,"adducthyp"))) {
      rdisop_res <- lapply(isomain, function(M0) {
        frag <- GetFragmentData(M0=M0, spec=spec, n=3, ionization=ionization)
        if (M0==isomain[length(isomain)]) {
          #browser()
          frag[1,] <- frag[1,]-attr(isomain,"adducthyp")+1.0073
        }
        rdisop_res <- EvaluateFragment(frag=frag, em=em, dppm=2*dppm, allowed_elements=allowed_elements, silent=silent, ionization=ionization, maxElements=maxElements, quick_isos=quick_isos, sf_db=sf_db)
        invisible(rdisop_res)
      })
    } else {
      rdisop_res <- lapply(isomain, function(M0) {
        frag <- GetFragmentData(M0=M0, spec=spec, n=3, ionization=ionization)
        rdisop_res <- EvaluateFragment(frag=frag, em=em, dppm=2*dppm, allowed_elements=allowed_elements, silent=silent, ionization=ionization, maxElements=maxElements, quick_isos=quick_isos, sf_db=sf_db)
        invisible(rdisop_res)
      })
    }
    RemoveEmptyFragments(rdisop_res, silent=silent, step="GetRdisopResult")
  }
  RemoveByScore <- function(rdisop_res, score_cutoff=0, silent=TRUE) {
    rdisop_res <- lapply(rdisop_res, function(x) { x[x[,"Score"]>=(score_cutoff*max(x[,"Score"])),] })
    rdisop_res <- RemoveEmptyFragments(rdisop_res, silent=silent, step="RemoveByScore")
  }
  GenerateMainOutput <- function(rdisop_res_list, stats, met_db) {
    cat(paste("\n\nTotal number of formulas per fragment before and after filtering...\n"))
    print(stats)
    rdisop_res_best <- rdisop_res_list[[1]]
    rdisop_res_best <- rdisop_res_best[!is.na(rdisop_res_best[,1]),]
    cat(paste("\n\nDetails of best candidate...\n"))
    print(rdisop_res_best)
    MH <- ifelse(nrow(rdisop_res_best)==1, 1, ifelse(abs(rdisop_res_best[nrow(rdisop_res_best),"Mass"]-rdisop_res_best[nrow(rdisop_res_best)-1,"Mass"]-72.0395)<0.0005, nrow(rdisop_res_best)-1, nrow(rdisop_res_best)))
    # check resulting candidate list against a DB for best matche
    if (!is.null(met_db)) {
      met_db[,"Formula"] <- enviPat::check_chemform(isotopes=as.matrix(allowed_elements,ncol=1), met_db[,"Formula"])[,"new_formula"]
      #M0 <- Rdisop::subMolecules(rdisop_res_best[MH,"Formula"], "H")$formula
      M0 <- enviPat::subform(rdisop_res_best[MH,"Formula"], "H1")
      if (any(met_db$Formula %in% M0)) {
        if (!silent) print(met_db[met_db$Formula %in% M0,])
        best_cand <- paste(met_db[met_db$Formula %in% M0, "Name"], collapse="; ")
        best_cand_col <- 3
      } else {
        if (any(abs(met_db$"M+H"-rdisop_res_best[MH,"Mass"])<0.002)) {
          if (!silent) print(met_db[abs(met_db$"M+H"-rdisop_res_best[MH,"Mass"])<0.002,-1])
          best_cand <- paste(met_db[abs(met_db$"M+H"-rdisop_res_best[MH,"Mass"])<0.002,"Name"],collapse="; ")
          best_cand_col <- 6
        } else {
          best_cand <- ""
          best_cand_col <- grDevices::grey(0.9)
        }
      }
    } else {
      best_cand <- ""
      best_cand_col <- 1
    }
    # outcommented to allow combination of spectrum plot with other plots
    # opar <- par(no.readonly=TRUE)
    # set colors for spectra plot (isomain peaks get a red color)
    tmp.col <- rep(1,nrow(spec))
    tmp.col[which(spec[,1] %in% isomain)] <- 2
    PlotSpec(x=spec, cols=tmp.col, txt=data.frame("mz"=rdisop_res_best[,"Mass"],"Formula"=rdisop_res_best[,"Formula"]), neutral_losses=neutral_losses, neutral_loss_cutoff=neutral_loss_cutoff, substitutions=substitutions, ionization=ionization)
    graphics::mtext(paste("Remaining combinations:",length(rdisop_res_list)), line=-1.2, adj=0, side=3, col=grDevices::grey(0.5))
    graphics::mtext(best_cand, line=-2.4, adj=0, side=3, col=best_cand_col)
    if (!is.null(correct_peak)) {
      graphics::mtext(correct_peak, line=-3.6, adj=0, side=3, col=grDevices::grey(0.5))
      fcor <- strsplit(correct_peak,", ")[[1]]
      fcor <- fcor[-1][grep("^C[[:digit:]]",fcor[-1])]
      if (length(fcor)==1) {
        fcor <- enviPat::check_chemform(isotopes=as.matrix(allowed_elements,ncol=1), chemforms=fcor)
        if (fcor[,1]) {
          warning("[InterpretMSSpectrum] Probably a wrong specification of 'correct_peak'", call. = FALSE)
        } else {
          #fcor <- Rdisop::addMolecules(fcor[,2], "H")$formula
          fcor <- enviPat::mergeform(fcor[,2], "H1")
        }
        # [Modification:if correct peak is specified but smaller than the 'observed' M0 the now uncommented line leads to trouble]
        # fcor <- which(sapply(rdisop_res_list, function(x) {fcor %in% x[,"Formula"]}))
        fcor <- which(sapply(rdisop_res_list, function(x) {fcor == x[nrow(x),"Formula"]}))
        cat(paste("\n\nRank of specified M0 =", ifelse(length(fcor)==1,fcor,"NA"), "\n"))
        graphics::mtext(paste0("Rank = ", fcor), line=-4.6, adj=0, side=3, col=grDevices::grey(0.5))
      }
    }
    # restore parameters
    # par(opar)
  }
  
  # read data (mz,int -table) from clipboard if not provided explicitly and sort by mz (should be standard but you never know...)
  if (is.null(spec)) spec <- ReadSpecClipboard()
  spec <- spec[order(spec[,1]),,drop=FALSE]
  
  # check  if spectrum is okay
  if (!nrow(spec)>=1) {
    
    #[ToDo] FurtherChecks may be usefull
    warning("[InterpretMSSpectrum] Spectra does not contain any information (nrow=0).", call. = FALSE)
    invisible(NULL)
    
  } else {
    
    # keep global error message for testing purposes (if correct_peak is known/provided)
    if (!is.null(correct_peak)) {
      #global_err <<- NULL # global assignment removed to pass RCheck without NOTEs
      global_err <- NULL
      local_check <- 0
    }
    #browser()
    
    # evaluate main peaks
    #browser()
    isomain <- DetermineIsomainPeaks(spec=spec, int_cutoff=0.03, precursor=precursor, ionization=ionization, ionmode=ionmode, limit=max_isomain_peaks)
    
    # # modify spectra if no good [M+H] is found but alternative adduct hypothesis instead
    # if (!is.null(attr(isomain,"adducthyp"))) {
    #   
    # }

    stats <- data.frame("mz"=round(isomain,4), "initial"=NA, "score_cutoff"=NA, "PlausibleFormula"=NA, "TypicalLosses"=NA)

    # start timing for testing purposes
    time_elapse <- Sys.time()
    # use Rdisop to get potential formulas (using up to n=2 isotopic peaks found)
    rdisop_res <- GetRdisopResult(spec=spec, isomain=isomain, silent=silent, em=em, dppm=dppm, allowed_elements=allowed_elements, ionization=ionization, maxElements=maxElements, quick_isos=quick_isos, sf_db=sf_db)
    stats[stats[,"mz"] %in% round(sapply(rdisop_res,attr,"M0"),4),"initial"] <- sapply(rdisop_res,nrow)
    if (!is.null(correct_peak)) {
      if (local_check==0 && length(grep(fml, rdisop_res[[length(rdisop_res)]][,1]))!=1) local_check <- 1
    }
    time_elapse <- c(time_elapse, Sys.time())
    
    # remove according to individual score based on mz deviation and isotopic fit
    rdisop_res <- RemoveByScore(rdisop_res, score_cutoff=score_cutoff, silent=silent)
    stats[stats[,"mz"] %in% round(sapply(rdisop_res,attr,"M0"),4),3] <- sapply(rdisop_res,nrow)
    
    if (!is.null(correct_peak) && local_check==0 && length(grep(fml, rdisop_res[[length(rdisop_res)]][,1]))!=1) local_check <- 2
    time_elapse <- c(time_elapse, Sys.time())
    
    # restrict to plausible formulas
    if (!sf_db) {
      # if a predetermined database is used there is no need to check formula plausibility again
      rdisop_res <- lapply(rdisop_res, function(x) { x[sapply(x[,"Formula"], PlausibleFormula, ruleset=ionization),] })
      rdisop_res <- RemoveEmptyFragments(rdisop_res, silent=silent, step="PlausibleFormula")
      stats[stats[,"mz"] %in% round(sapply(rdisop_res,attr,"M0"),4),4] <- sapply(rdisop_res,nrow)
    } else {
      #stats[stats[,"mz"] %in% round(sapply(rdisop_res,attr,"M0"),4),4] <- sapply(rdisop_res,nrow)
    }
    
    if (!is.null(correct_peak) && local_check==0 && length(grep(fml, rdisop_res[[length(rdisop_res)]][,1]))!=1) local_check <- 3
    time_elapse <- c(time_elapse, Sys.time())
    #rdisop_res <<- rdisop_res
    
    # restrict based on neutral losses (only if more than 1 fragment is left in spectrum)
    # test of all losses are potentially helpfull
    nl_vec <- sapply(neutral_losses[,"Formula"], function(x) { neutral_losses[neutral_losses[,"Formula"]==x,"Mass"] })
    if (length(rdisop_res)>=2) {
      for (k in 1:length(nl_vec)) {
        rdisop_res <- RestrictByTypicalLosses(rdisop_res=rdisop_res, tl=nl_vec[k], neutral_loss_cutoff=neutral_loss_cutoff, punish=0.5, substitutions=substitutions)
        if (!is.null(correct_peak) && local_check==0 && length(grep(fml, rdisop_res[[length(rdisop_res)]][,1]))!=1) {
          local_check <- 4
          #print(paste(names(nl_vec)[k], "caused a problem."))
        }
      }
    }
    stats[stats[,"mz"] %in% round(sapply(rdisop_res,attr,"M0"),4),5] <- sapply(rdisop_res,nrow)
    time_elapse <- c(time_elapse, Sys.time())

    # close SQLite connection
    if (sf_db) DBI::dbDisconnect(db_con)
    
    if (sum(sapply(rdisop_res,nrow))>0) {
      
      # obtain most likely combination
      #rdisop_res <<- rdisop_res
      #browser()
      rdisop_res_list <- ScoreFormulaCombination(rdisop_res, nl_vec=nl_vec, punish_invalid=0.5, punish_S=0.2, punish_nonplausible=0.5, return_rank=NA, neutral_loss_cutoff=neutral_loss_cutoff, substitutions=substitutions, silent=silent)
      time_elapse <- c(time_elapse, Sys.time())
      
      # plot annotated spectrum
      if (!silent) GenerateMainOutput(rdisop_res_list, stats, met_db)
      time_elapse <- c(time_elapse, Sys.time())
      
      # write potential error source to global
      if (!is.null(correct_peak) && local_check>0) {
        global_err <- paste(correct_peak, "not present after", c("initial generation","score filter","plausibility filter","neutral loss filter")[local_check])
        #assign("global_err", global_err, envir=.GlobalEnv)
        print(global_err)
      }
      
      time_elapse <- round(diff(time_elapse),4); names(time_elapse) <- c("Rdisop","ScoreFilt","Plausible","NeutralLoss","PathEval","Plot")
      if (!silent) {
        cat("\n\nTime elapsed during individual processing steps...\n")
        print(time_elapse)
      }
      
      # return final result list
      attr(rdisop_res_list,"stats") <- stats
      invisible(rdisop_res_list)
      
    } else {
      
      invisible(NULL)
      
    }
  }
}
