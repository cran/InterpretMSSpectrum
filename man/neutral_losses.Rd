\name{neutral_losses}
\alias{neutral_losses}
\docType{data}
\title{
neutral_losses
}
\description{
A data table defining neutral losses in GC-APCI-MS for silylated compounds.
}
\usage{data("neutral_losses")}
\format{
  A data frame with 22 observations on the following 3 variables.
  \describe{
    \item{\code{Name}}{a character vector}
    \item{\code{Formula}}{a character vector}
    \item{\code{Mass}}{a numeric vector}
  }
}
\details{
The data frame consists 2 character columns ('Name' and 'Formula') and the numeric column 'Mass'.
In a mass spectrum peak pairs are analyzed for mass differences similar to the ones defined in neutral_losses.
If such a mass difference is observed we can assume that the according 'Formula' is the true neutral loss
observed in this spectrum. In a plot this peak pair would be connected by a grey line and annotated with
the information from 'Name'.
}
\source{
This list has been put together mnually by Jan Lisec analyzing several GC-APCI-MS data sets.
}
\references{
See published paper in $$
}
\examples{
data(neutral_losses)
str(neutral_losses)
}
\keyword{datasets}
