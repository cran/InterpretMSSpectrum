#' @name sendToMSF
#' @aliases sendToMSF
#' @aliases sendToMSF.default
#' @aliases sendToMSF.findMAIN
#' 
#' @title Exporting spectra to MSFinder.
#' 
#' @description Send spectrum to MSFinder.
#' 
#' @details 
#' In the default case 'x' can be a matrix or data frame, where the first two columns are assumed to contain the 'mz' and 'intensity' values, respectively. Further arguments 'precursormz' and 'precursortype' are required in this case. Otherwise 'x' can be of class \code{findMAIN}.
#' 
#' @param x A matrix or 'findMAIN' object
#' @param ... Arguments passed to methods of \code{\link{writeMSF}}.
#'
#' @return
#' Full path of generated MAT file (invisibly).
#' 
#' @rdname sendToMSF
#' @export
#' 
#' @references H.Tsugawa et al (2016) Hydrogen rearrangement rules: computational MS/MS fragmentation and structure elucidation using MS-FINDER software. Analytical Chemistry, 88, 7946-7958
#' 
sendToMSF <- function(x, ...) {
  UseMethod("sendToMSF")
}

#' @param precursormz m/z of (de)protonated molecule or adduct ion 
#' @param precursortype adduct type, e.g. "[M+H]+", "[M+Na]+". Accepted values are all adduct ions supported by MSFINDER.
#' @param outfile Name of MAT file. If NULL, a temporary file is created in the  per-session temporary directory (see \code{\link{tempdir}}).
#' @param MSFexe Full path of MS-FINDER executable. This needs to be set  according to your system. If \code{NULL}, MAT files are written but the program is not opened.
#'   
#' @rdname sendToMSF
#' @method sendToMSF default
#' 
sendToMSF.default <- function(x, 
                              precursormz,
                              precursortype="[M+H]+",
                              outfile=NULL, 
                              MSFexe=NULL, 
                              ...) {
  if(is.null(outfile)) {
    fn <- sprintf("unknown%s_%.2f.mat", 
                  format(Sys.time(), "%Y%m%d_%H%M%OS2"),
                  x[,1][which.max(x[,2])])
    outfile <- file.path(tempdir(), fn) # use the per-session temporary dir
  }
  writeMSF.default(x, 
                   precursormz=precursormz, 
                   precursortype=precursortype, 
                   outfile=outfile, ...)
  if(!is.null(MSFexe))
    system2(normalizePath(MSFexe), normalizePath(dirname(outfile)), wait=FALSE)
  invisible(outfile)
}

#' @rdname sendToMSF
#' @method sendToMSF findMAIN
#' 
#' @param rank Which rank from 'findMAIN' should be exported
#' @param ms2spec An (optional) MS2 spectrum to be passed to MSFINDER. If  \code{NULL}, the MS1 spectrum used by 'findMAIN' is used. If dedicated MS2 spectra are available, this option should be used.
#'   
sendToMSF.findMAIN <- function(x, 
                               rank=1,
                               ms2spec=NULL,
                               outfile=NULL, 
                               MSFexe=NULL, 
                               ...) {
  if(is.null(outfile)) {
    spec <- x[[rank]][,1:2,drop=F]
    fn <- sprintf("unknown%s_%.2f.mat", 
                  format(Sys.time(), "%Y%m%d_%H%M%OS2"),
                  attr(x[[rank]], "scores")[, "prec"])
    outfile <- file.path(tempdir(), fn) # use the per-session temporary dir
  }
  writeMSF.findMAIN(fmo=x, 
                    rank=rank,
                    ms2spec=ms2spec,
                    outfile=outfile, 
                    ...)
  if(!is.null(MSFexe)) system2(normalizePath(MSFexe), normalizePath(dirname(outfile)), wait=FALSE)
  invisible(outfile)
}