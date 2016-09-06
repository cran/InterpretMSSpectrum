\name{test_spectrum}
\alias{test_spectrum}
\docType{data}
\title{
test_spectrum
}
\description{
Example spectrum of Glutamic acid (3TMS) measured on a Bruker impact II.
}
\usage{data("test_spectrum")}
\format{
  A data frame with 47 observations on the following 2 variables.
  \describe{
    \item{\code{mz}}{a numeric vector}
    \item{\code{int}}{a numeric vector}
  }
}
\examples{
data(test_spectrum)
head(test_spectrum)
PlotSpec(test_spectrum)
}
\keyword{datasets}
\keyword{internal}