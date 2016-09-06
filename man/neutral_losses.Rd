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
This list has been put together mnually by Jan Lisec analyzing multiple GC-APCI-MS data sets.
}
\examples{
data(neutral_losses)
str(neutral_losses)
}
\keyword{datasets}
\keyword{internal}