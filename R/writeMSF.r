#' @name writeMSF
#' @aliases writeMSF
#' @aliases writeMSF.default
#' @aliases writeMSF.findMAIN
#' 
#' @title writeMSF.
#' 
#' @description Write a spectrum file in MSFinder format.
#' 
#' @details 
#' In the default case 'x' can be a matrix or data frame, where the first two columns are assumed to contain the 'mz' and 'intensity' values, respectively. Further arguments 'precursormz' and 'precursortype' are required in this case. Otherwise 'x' can be of class \code{findMAIN}.
#' 
#' @param x A matrix, 'findMAIN' or other object for which methods are defined.
#' @param ... Arguments passed to method writeMSF.findMAIN.
#'
#' @return
#' Write spectrum to MAT file for evaluation in MSFinder
#' 
#' @rdname writeMSF
#' @export
#' @keywords internal
#' 
writeMSF <- function(x, ...) {
    UseMethod("writeMSF")
}

#' @param precursormz Precursor m/z
#' @param precursortype Precursor type
#' @param name Give the spectrum a name
#' @param ionmode "Positive" or "Negative"
#' @param ms1spec (Optional) MS1 spectrum. 
#' @param retentiontime (Optional) retention time of the spectrum that will be used by MSFinder for refined prediction.
#' @param outfile Name of MAT file, or \code{NULL} for 'stdout'.
#'
#' @rdname writeMSF
#' @export
#'
writeMSF.default <- function(x, 
                             precursormz, 
                             precursortype="[M+H]+", 
                             name="unknown",
                             ionmode="Positive", 
                             ms1spec=NULL, 
                             retentiontime=NULL, 
                             outfile=NULL,
                             ...) {
  "%&%" <- function(x, y) paste0(x, y)
  out <- "NAME: " %&% name %&% "\n"
  if (!is.null(retentiontime)) {
    out <- out %&% "RETENTIONTIME: " %&% retentiontime %&% "\n"
  }
  out <- out %&% "PRECURSORMZ: " %&% sprintf("%.4f", precursormz) %&% "\n"
  out <- out %&% "PRECURSORTYPE: " %&% precursortype %&% "\n"
  out <- out %&% "IONMODE: " %&% ionmode %&% "\n"
  if (!is.null(ms1spec) && nrow(ms1spec)>0 && any(ms1spec[,2]>0) && any(ms1spec[,2]>0)) {
    out <- out %&% "MSTYPE: MS1\n"
    out <- out %&% "Num Peaks: " %&% nrow(ms1spec[ms1spec[,2]>0,,drop=F]) %&% "\n"
    for (i in which(ms1spec[,2]>0)) out <- out %&% ms1spec[i, 1] %&% " " %&% ms1spec[i, 2] %&% "\n"
  }
  if (nrow(x)>0 && any(!is.na(x[,2])) && any(x[,2]>0)) {
    out <- out %&% "MSTYPE: MS2\n"
    out <- out %&% "Num Peaks: " %&% nrow(x[x[,2]>0,,drop=F]) %&% "\n"
    for (i in which(x[,2]>0)) out <- out %&% x[i, 1] %&% " " %&% x[i, 2] %&% "\n"
  }
  cat(out, file=if(is.null(outfile)) "" else outfile)
}

#' @param rank Which rank from 'findMAIN' should be exported
#' @param ms2spec If available you may provide the according MS2 spectrum
#'
#' @rdname writeMSF
#' @export
#'
writeMSF.findMAIN <- function(x, 
                              rank=1, 
                              ms2spec=NULL, 
                              outfile=NULL, 
                              ...) {
  precursors_tested <- attr(x, "adducthyp_tested")
  precursordict <- data.frame(massdiff=unlist(sapply(precursors_tested, getRuleFromIonSymbol)[4,]),
                              name=precursors_tested,
                              stringsAsFactors=FALSE)
  nranks <- length(x)
  r <- rank
  r[r>nranks] <- nranks
  r <- unique(r)
  for(rnk in r) {
    ms1spec <- x[[1]][,1:2, drop=F]
    ms2spec <- if(is.null(ms2spec)) {
      ms1spec 
    } else if(is.list(ms2spec) && !is.data.frame(ms2spec)) {
      ms2spec[[rnk]] 
    } else {
      ms2spec 
    }
    scoresObj <- attr(x[[rnk]], "scores")
    precursormz <- scoresObj[,"adductmz"]
    adducthyp <- scoresObj[,"adducthyp"]
    precursortype <- precursordict[,2][which(abs(adducthyp-precursordict[,1])<0.01)]
    writeMSF.default(x=ms2spec, 
                     precursormz=precursormz, 
                     precursortype=precursortype,
                     ms1spec=ms1spec,
                     name=if(is.null(outfile)) "unknown" else {
                       if(length(r)==1) outfile else sprintf("%s_%02d", sub("\\.mat", "", basename(outfile)), rnk)
                     },
                     outfile=if(is.null(outfile)) "" else {
                       if(length(r)==1) outfile else file.path(dirname(outfile), sprintf("%s_%02d.mat", sub("\\.mat", "", basename(outfile)), rnk))
                     },
                     ...)
  }
}
