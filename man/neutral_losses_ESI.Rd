\name{neutral_losses_ESI}
\alias{neutral_losses_ESI}
\docType{data}
\title{
neutral_losses_ESI
}
\description{
A data table defining neutral losses in LC-ESI-MS (positive mode).
}
\usage{data("neutral_losses_ESI")}
\format{
  A data frame with 45 observations on the following 3 variables.
  \describe{
    \item{\code{Name}}{a character vector containing the fragment name used for plot annnotation}
    \item{\code{Formula}}{a character vector containing chemical formulas}
    \item{\code{Mass}}{a numeric vector containing the mass according to Formula}
  }
}
\details{
The data frame consists of two character columns ('Name' and 'Formula') and the numeric column 'Mass'.
In a mass spectrum peak pairs are analyzed for mass differences similar to the ones defined in neutral_losses.
If such a mass difference is observed, we can assume that the according 'Formula' is the true neutral loss
observed in this spectrum. In a plot this peak pair would be connected by a grey line and annotated with
the information from 'Name'. In formula evaluation this peak pair would be used to limit formula suggestions
with respect to plausability, i.e. if mass fragments A and B exist with mass difference 16.0313 than we can
assume that the respective sum formulas have to be different by CH4. In consequence we can exlude sum formula
suggestions for B which do not have a valid corresponding sum formula in A and vice versa.
}
\source{
This list has been put together mnually by Jan Lisec analyzing multiple LC-ESI-MS (positive mode) data sets.
}
\examples{
data(neutral_losses_ESI)
str(neutral_losses_ESI)
}
\keyword{datasets}
\keyword{internal}