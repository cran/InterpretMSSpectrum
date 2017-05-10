#'@title IMS_parallel.
#'
#'@description
#'\code{IMS_parallel} is a parallel implementation of \code{\link{InterpretMSSpectrum}}.
#'
#'@details
#'For mass processing and testing it may be sufficient to use InterpretMSSpectrum without plotting functionality. However, function is likely to be deprecated or integrated as an option into the main function in the future.
#'
#'@param spectra List of spectra.
#'@param precursor vector of precursor masses of length(spectra).
#'@param ncores Number of cores available.
#'@param correct_peak Potentially a vector of correct Peaks, see \code{InterpretMSSpectrum} for details.
#'@param ... Further parameters passed directly to \code{InterpretMSSpectrum}.
#'
#'@return
#'A list of \code{InterpretMSSpectrum} result objects which can be systematically evaluated. However, note that plotting is unfortunately not enabled for parallel processing.
#'
#'@seealso\code{\link{InterpretMSSpectrum}}
#'
#'@import foreach
#'@import doParallel
#'
#'@export
#'
IMS_parallel <- function(spectra=NULL, ncores=8, precursor=NULL, correct_peak=NULL, ...) {
  # dummy line to initialize idx
  idx <- NULL
  cat(paste0("\n\nStart processing of ", length(spectra), " spectra, which may take a while.\n..."))
  subset <- split(1:length(spectra), gl(ncores,ceiling(length(spectra)/ncores)))
  cl <- parallel::makeCluster(ncores)
  doParallel::registerDoParallel(cl)
  res <- foreach::foreach(idx = subset)  %dopar%  {
    lapply(idx, function(i) {
      cat(correct_peak[i])
      InterpretMSSpectrum::InterpretMSSpectrum(spec=spectra[[i]], silent=TRUE, precursor=precursor[i], correct_peak=correct_peak[i], ...)
    })
  }
  doParallel::stopImplicitCluster()
  invisible(unlist(res, recursive = FALSE, use.names=FALSE))
}