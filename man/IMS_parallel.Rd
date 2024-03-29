% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IMS_parallel.R
\name{IMS_parallel}
\alias{IMS_parallel}
\title{IMS_parallel.}
\usage{
IMS_parallel(
  spectra = NULL,
  ncores = 8,
  precursor = NULL,
  correct_peak = NULL,
  ...
)
}
\arguments{
\item{spectra}{List of spectra.}

\item{ncores}{Number of cores available.}

\item{precursor}{vector of precursor masses of length(spectra).}

\item{correct_peak}{Potentially a vector of correct Peaks, see \code{InterpretMSSpectrum} for details.}

\item{...}{Further parameters passed directly to \code{InterpretMSSpectrum}.}
}
\value{
A list of \code{InterpretMSSpectrum} result objects which can be systematically evaluated. 
    However, note that plotting is unfortunately not enabled for parallel processing.
}
\description{
\code{IMS_parallel} is a parallel implementation of \code{\link{InterpretMSSpectrum}}.
}
\details{
For mass processing and testing it may be sufficient to use \code{InterpretMSSpectrum} 
    without plotting functionality. However, function is likely to be deprecated or integrated 
    as an option into the main function in the future.
}
\seealso{
\code{\link{InterpretMSSpectrum}}
}
