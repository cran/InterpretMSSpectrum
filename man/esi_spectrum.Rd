\name{esi_spectrum}
\alias{esi_spectrum}
\docType{data}
\title{
esi_spectrum
}
\description{
Example spectrum of Glutamic acid (3TMS) measured on a Bruker impact II.
}
\usage{data("esi_spectrum")}
\format{
  A matrix with 42 observations on the following 2 variables.
  \describe{
    \item{\code{mz}}{a numeric vector}
    \item{\code{int}}{a numeric vector}
  }
}
\examples{
data(esi_spectrum)
head(esi_spectrum)
PlotSpec(esi_spectrum)
}
\keyword{datasets}
\keyword{internal}