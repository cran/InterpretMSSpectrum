% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getRuleFromIonSymbol.R
\name{getRuleFromIonSymbol}
\alias{getRuleFromIonSymbol}
\title{Generate adduct rule from ion symbol}
\usage{
getRuleFromIonSymbol(ions = "[M+H]+")
}
\arguments{
\item{ions}{character vector of ion symbols, e.g. "[M+H]+", "[M+Na]+", 
"[M+H-NH3]-". Please use full notation in square brackets, though some 
frequent ions can be abbrevated ("M+H","M+Na","M+K","M+NH4", "M+", "M", 
"M-H","M+Cl-", "M-").}
}
\value{
A data frame with four columns "name", "nmol", "charge", "massdiff".
}
\description{
Translate an ion symbol to an adduct rule. This function is used internally by \code{findMAIN}, but may be useful elsewhere.
}
\examples{
getRuleFromIonSymbol(c("[M+H]+", "[M+Na]+"))
}
\keyword{internal}
