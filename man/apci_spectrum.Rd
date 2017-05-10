\name{apci_spectrum}
\alias{apci_spectrum}
\docType{data}
\title{
apci_spectrum
}
\description{
Example spectrum of Glutamic acid (3TMS) measured on a Bruker impact II.
}
\usage{data("apci_spectrum")}
\format{
  A data frame with 47 observations on the following 2 variables.
  \describe{
    \item{\code{mz}}{a numeric vector}
    \item{\code{int}}{a numeric vector}
  }
}
\examples{
data(apci_spectrum)
head(apci_spectrum)
PlotSpec(apci_spectrum)
}
\keyword{datasets}
\keyword{internal}