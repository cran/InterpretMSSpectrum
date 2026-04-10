#'@title Plot Mass Spectrum.
#'
#'@description
#'\code{PlotSpec} will read, evaluate and plot a deconvoluted mass spectrum (mass*intensity pairs) from TMS-derivatized GC-APCI-MS data.
#'The main purpose is to visualize the relation between deconvoluted masses.
#'
#'@param x A two-column matrix with ("mz", "int") information.
#'@param masslab The cutoff value  (relative to basepeak) for text annotation of peaks.
#'@param rellab TRUE/FALSE. Label masses relative to largest mass in plot (if TRUE), 
#'    absolute (if FALSE) or to specified mass (if numeric).
#'@param cutoff Show only peaks with intensity higher than cutoff*I(base peak). 
#'    This will limit the x-axis accordingly.
#'@param cols Color vector for peaks with length(cols)==nrow(x).
#'@param txt Label peaks with specified text (column 1 specifies x-axis value, 
#'    column 2 specifies label). Labels will be converted to expressions following
#'    the ruleset to format chemical formulas (see Examples).
#'@param mz_prec Numeric precision of m/z (=number of digits to plot).
#'@param neutral_losses Data frame of defined building blocks (Name, Formula, Mass). 
#'    If not provided data("neutral_losses") will be used.
#'@param neutral_loss_cutoff Specifies the allowed deviation in mDa for neutral 
#'    losses to be accepted from the provided neutral loss list.
#'@param substitutions May provide a two column table of potential substitutions 
#'    (for adducts in ESI-MS).
#'@param ionization Either APCI or ESI (important for main peak determination).
#'@param precursor Internally main peaks will be determined up to a supposed 
#'    precursor obtained by `DetermineIsomainPeaks` and annotations will only be 
#'    plotted up to this mass. To plot annotations for the full mass range, set 
#'    `precursor` to a higher mass.
#'@param xlim To specify xlim explicitly (for comparative plotting).
#'@param ylim To specify ylim explicitly (for comparative plotting).
#'@param mar Set to NULL to allow specifying margins outside the function (important for use in Shiny App context with brush event).
#'
#'@return
#'An annotated plot of the mass spectrum.
#'
#'@examples
#'# load APCI test data and apply function
#'utils::data(apci_spectrum, package = "InterpretMSSpectrum")
#'PlotSpec(x=apci_spectrum, ionization="APCI")
#'
#'# normalize test data by intensity
#'s <- apci_spectrum
#'s[,2] <- s[,2]/max(s[,2])
#'PlotSpec(x=s)
#'
#'# use relative labelling
#'PlotSpec(x=s, rellab=s[which.max(s[,2]),1])
#'
#'# avoid annotation of masses and fragments
#'PlotSpec(x=s, masslab=NULL, neutral_losses=NA)
#'
#'# provide individual neutral loss set
#'tmp <- data.frame("Name"=c("Loss1","C2H6"),"Formula"=c("",""),"Mass"=c(90.05,27.995))
#'PlotSpec(x=s, neutral_losses=tmp)
#'
#'# provide additional color and annotation information per peak
#'PlotSpec(x=s, cols=1+(s[,2]>0.1), txt=data.frame("x"=s[s[,2]>0.1,1],"txt"="txt"))
#'
#'# annotate a sum formula
#'PlotSpec(x=s, txt=data.frame("x"=s[which.max(s[,2]),1],"txt"="C6H12O6"))
#'
#'# or annotate a mix of sum formula and text
#'txt <- data.frame(
#'  "x" = s[order(s[,2],decreasing=TRUE)[1:3],1],
#'  "txt" = c("C6H12O6", "some text", "H2SO3")
#')
#'PlotSpec(x=s, txt=txt)
#'
#'# if 'txt' contains a column SMILES with valid SMILES codes, they will be rendered additionally
#'txt[,"SMILES"] <- c("O=CC(O)C(O)C(O)C(O)CO", "", "OS(=O)(=O)")
#'PlotSpec(x=s, txt=txt)
#'txt[,"w"] <- rep(0.4,3)
#'PlotSpec(x=s, txt=txt)
#'
#'# simulate a Sodium adduct to the spectrum (and annotate using substitutions)
#'p <- which.max(s[,2])
#'s <- rbind(s, c(21.98194+s[p,1], 0.6*s[p,2]))
#'PlotSpec(x=s, substitutions=matrix(c("H","Na"),ncol=2,byrow=TRUE))
#'
#'# add structural formulas by providing SMILES codes
#'
#'# use predefined ESI parameter set
#'PlotSpec(x=s, ionization="ESI")

#'@export

PlotSpec <- function(x=NULL, masslab=0.1, rellab=FALSE, cutoff=0.01, cols=NULL, txt=NULL, mz_prec=4, ionization=NULL, neutral_losses=NULL, neutral_loss_cutoff=NULL, substitutions=NULL, precursor=NULL, xlim=NULL, ylim=NULL, mar=c(2,2,0.5,0)+0.5) {
  # potential parameters
  max_isomain_peaks <- NULL

  #
  if (is.null(neutral_loss_cutoff)) {
    if (is.null(ionization)) {
      neutral_loss_cutoff <- 1
    } else {
      if (ionization=="APCI") neutral_loss_cutoff <- 0.5
      if (ionization=="ESI") neutral_loss_cutoff <- 2
    }
  }

  # check spectra format
  x <- x[,1:2,drop=FALSE]
  test_valid_spectra <- all(apply(x,2,is.numeric)) & prod(dim(x))>1

  if (test_valid_spectra) {

    # load building block definition
    if (is.null(neutral_losses)) {
      if (is.null(ionization)) {
        neutral_losses <- NULL
      } else {
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
      }
    }

    # attach substitution information to neutral loss table
    if (!is.null(substitutions)) {
      substitutions <- data.frame(
        "Name" = apply(substitutions, 1, paste, collapse="/"), 
        "Formula" = rep("", nrow(substitutions)), 
        "Mass" = apply(substitutions, 1, function(x) {
          abs(diff(get_exactmass(x)))
        })
      )
      if (!exists("neutral_losses", envir = environment())) {
        neutral_losses <- substitutions
      } else {
        neutral_losses <- rbind(neutral_losses[, colnames(substitutions)], substitutions)
      }
      # remove duplicated entries
      neutral_losses <- neutral_losses[!duplicated(neutral_losses[,"Mass"]),,drop=FALSE]
    }

    # define mass range for plotting
    xf <- rep(T, length(x[,1]))
    if (is.numeric(cutoff)) xf[x[,2] < cutoff*max(x[,2])] <- FALSE

    # cols will be provided by InterpretMSSpec
    if (is.null(cols)) cols <- rep(1, nrow(x))

    # set up main plot
    if (!is.null(mar) && is.numeric(mar) && length(mar)==4L && all(is.finite(mar))) {
      opar <- graphics::par("mar")
      on.exit(graphics::par("mar" = opar))
      graphics::par("mar"=mar)
    }
    if (is.null(xlim)) {
      xlim <- c(floor(min(x[xf,1])), ceiling(max(x[xf,1])))
      if (diff(xlim)<30) xlim <- xlim + c(-1,1)*ceiling((30-diff(xlim))/2)
    }
    #xlim <- round(xlim)
    if (is.null(ylim)) {
      ylim <- range(c(0, x[xf,2]))
    }
    graphics::plot(x=x[xf,1], y=x[xf,2], type="h", las=1, xlab="", ylab="", main="", col=cols[xf], ann=F, axes=F, xlim=xlim, ylim=ylim)
      graphics::axis(2)
      graphics::axis(1, tcl=-0.8, lwd=1.2)
      if (diff(xlim)<600) graphics::axis(1, at=seq(floor(min(xlim)/10)*10,ceiling(max(xlim)/10)*10,10), labels=FALSE, tcl=-0.6, lwd=1.2)
      if (diff(xlim)<100) graphics::axis(1, at=seq(min(xlim),max(xlim),1), labels=FALSE, tcl=-0.3)
      graphics::box()

    # indicate typical losses
    if (prod(dim(neutral_losses))>1) {
      # determine the main peaks of all isotopic clusters
      isomain <- which(x[xf,1] %in% DetermineIsomainPeaks(spec=x[xf,1:2,drop=F], int_cutoff=0.03, ionization=ifelse(is.null(ionization),"APCI",ionization), limit=max_isomain_peaks, precursor = precursor))
      # get distance matrix for isomain peaks and annotate typical losses
      dmz <- sapply(x[xf,1][isomain], function(y) {y-x[xf,1][isomain]})
      for (i in 1:nrow(neutral_losses)) {
        l <- which(abs(dmz[upper.tri(dmz)]-neutral_losses[i,3])<=(neutral_loss_cutoff/1000))
        if (length(l)>0) {
          l <- cbind(row(dmz)[upper.tri(dmz)][l],col(dmz)[upper.tri(dmz)][l])
          for (j in 1:nrow(l)) {
            graphics::lines(x[xf,][isomain,][l[j,],], lty=2, col=grDevices::grey(0.8))
            if (!inherits(try(fml2expr(neutral_losses[i,1]), silent=TRUE), "try-error")) {
              graphics::text(x=mean(x[xf,][isomain,][l[j,],1]), y=mean(x[xf,][isomain,][l[j,],2]), labels=fml2expr(neutral_losses[i,1]), col=grDevices::grey(0.8), cex=graphics::par("cex")*0.8)  
            } else {
              graphics::text(x=mean(x[xf,][isomain,][l[j,],1]), y=mean(x[xf,][isomain,][l[j,],2]), labels=neutral_losses[i,1], col=grDevices::grey(0.8), cex=graphics::par("cex")*0.8)
            }
          }
        }
      }
    }

    # print text (sum formulas and SMILES where provided)
    if (!is.null(txt)) {
      tmp.y <- sapply(txt[,1],function(y){x[which.min(abs(x[,1]-y)),2]})
      if (!inherits(try(fml2expr(txt[,2]), silent=TRUE), "try-error")) {
        for (i in 1:nrow(txt)) graphics::text(x=txt[i,1], y=tmp.y[i], pos=ifelse(tmp.y[i]>0.9*max(x[xf,2]),1,3), labels=fml2expr(txt[i,2]), col=4, cex=graphics::par("cex")*0.8)
      } else {
        graphics::text(x=txt[,1], y=tmp.y, pos=sapply(tmp.y, function(y) {ifelse(y>0.9*max(x[xf,2]),1,3)}), labels=txt[,2], col=4, cex=graphics::par("cex")*0.8)
      }
      if ("SMILES" %in% colnames(txt)) {
        gp_usr <- c(xlim+c(-1,1)*0.039*diff(xlim), ylim+c(-1,1)*0.039*diff(ylim))
        for (i in which(txt[,"SMILES"]!="")) {
          coords <- square_subplot_coord(x = txt[i,1], y = tmp.y[i], gp_usr = gp_usr, w = ifelse("w" %in% colnames(txt), txt[i,"w"], 0.25))
          renderSMILES(smiles = txt[i,"SMILES"], coords = coords, gp_usr = gp_usr)
        }
      }
    }

    # print masses
    if (is.numeric(masslab) && length(masslab)==1 && masslab>=0 && masslab<=1) {
      xf <- xf & (x[,2] >= masslab*max(x[xf,2]))
      if (length(rellab)==1 && is.numeric(rellab)) tmp.lab <- round(x[xf,1]-rellab, mz_prec)
      if (length(rellab)==1 && is.logical(rellab)) tmp.lab <- round(x[xf,1]-ifelse(rellab, x[xf,1][which.max(x[xf,2])], 0), mz_prec)
      if (exists("tmp.lab") && length(tmp.lab)==sum(xf)) graphics::text(x=x[xf,1], y=x[xf,2], labels=tmp.lab, cex=graphics::par("cex")*0.7)
    }

  } else {

    plot(1, 1, axes=F, ann=F, type="n")
    graphics::text(x=1, y=1, "No valid Spectra format...")

  }

  invisible(NULL)
}
