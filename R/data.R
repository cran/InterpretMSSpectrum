#' @title Default adduct lists used by \code{findMAIN}.
#' @docType data
#' @format A list of two character vectors:
#' \describe{
#' \item{Positive}{default adducts used in ESI(+) mode}
#' \item{Negative}{default adducts used in ESI(-) mode}
#' }
#' @usage data(Adducts)
#' @source A reasonable selection of frequent adducts based on the list in R-package CAMERA
"Adducts"

#' @title List of chemical elements.
#' @docType data
#' @format A data frame with 103 observations on the following 2 variables:
#' \describe{
#' \item{name}{a character vector of elemental abbreviations}
#' \item{mass}{a numeric vector of exact masses of the elements main isotope}
#' }
#' @usage data(chemical_elements)
"chemical_elements"

#' @title Orbitrap spectra
#' @docType data
#' @format A list with 550 matrices (spectra). Two attributes are attached to 
#' each spectrum: 
#' \describe{
#' \item{Formula}{sum formula of (neutral) compound}
#' \item{ExactMass}{exact mass of (neutral) compound}
#' }
#' @usage data(OrbiMS1)
#' @description A set of 550 MS1 pseudo-spectra of metabolite standards, acquired on an 
#' Orbitrap-type mass analyzer (Q Exactive, Thermo-Fisher) in electrospray ionization (ESI) 
#' positive mode. Spectra were generated from Thermo raw files using xcms/CAMERA.
"OrbiMS1"

#' @title APCI spectrum
#' @docType data
#' @format A data frame with 47 observations on the following 2 variables: 
#' \describe{
#' \item{mz}{a numeric vector ion masses}
#' \item{int}{a numeric vector of intensities}
#' }
#' @usage data(apci_spectrum)
#' @description Example spectrum of Glutamic acid (3TMS) measured on a Bruker impact II.
"apci_spectrum"

#' @title ESI spectrum
#' @docType data
#' @format A data frame with 42 observations on the following 2 variables: 
#' \describe{
#' \item{mz}{a numeric vector ion masses}
#' \item{int}{a numeric vector of intensities}
#' }
#' @usage data(esi_spectrum)
#' @description Example spectrum of Sedoheptulose 7-phosphate measured on a Bruker impact II.
"esi_spectrum"

#' @title A data table defining neutral losses in LC-ESI-MS (positive mode).
#' @docType data
#' @format A data frame with 45 observations on the following 3 variables:
#'  \describe{
#'    \item{\code{Name}}{a character vector containing the fragment name used for plot annnotation}
#'    \item{\code{Formula}}{a character vector containing chemical formulas}
#'    \item{\code{Mass}}{a numeric vector containing the mass according to Formula}
#'  }
#' @details The data frame consists of two character columns ('Name' and 'Formula') and the numeric column 'Mass'.
#' In a mass spectrum peak pairs are analyzed for mass differences similar to the ones defined in neutral_losses.
#' If such a mass difference is observed, we can assume that the according 'Formula' is the true neutral loss
#' observed in this spectrum. In a plot this peak pair would be connected by a grey line and annotated with
#' the information from 'Name'. In formula evaluation this peak pair would be used to limit formula suggestions
#' with respect to plausibility, i.e. if mass fragments A and B exist with mass difference 16.0313 than we can
#' assume that the respective sum formulas have to be different by CH4. In consequence we can exclude sum formula
#' suggestions for B which do not have a valid corresponding sum formula in A and vice versa.
#' @source This list has been put together manually by Jan Lisec analyzing multiple LC-ESI-MS (positive mode) data sets.
#' @usage data(neutral_losses_ESI)
#' @keywords datasets
"neutral_losses_ESI"

#' @title A data table defining typical neutral losses in GC-APCI-MS for silylated compounds.
#' @docType data
#' @format A data frame with 22 observations on the following 3 variables:
#'  \describe{
#'    \item{\code{Name}}{a character vector containing the fragment name used for plot annnotation}
#'    \item{\code{Formula}}{a character vector containing chemical formulas}
#'    \item{\code{Mass}}{a numeric vector containing the mass according to Formula}
#'  }
#' @details The data frame consists of two character columns ('Name' and 'Formula') and the numeric column 'Mass'.
#' In a mass spectrum peak pairs are analyzed for mass differences similar to the ones defined in neutral_losses.
#' If such a mass difference is observed, we can assume that the according 'Formula' is the true neutral loss
#' observed in this spectrum. In a plot this peak pair would be connected by a grey line and annotated with
#' the information from 'Name'. In formula evaluation this peak pair would be used to limit formula suggestions
#' with respect to plausibility, i.e. if mass fragments A and B exist with mass difference 16.0313 than we can
#' assume that the respective sum formulas have to be different by CH4. In consequence we can exclude sum formula
#' suggestions for B which do not have a valid corresponding sum formula in A and vice versa.
#' @source This list has been put together manually by Jan Lisec analyzing multiple GC-APCI-MS data sets.
#' @usage data(neutral_losses_ESI)
#' @keywords datasets
"neutral_losses_APCI"

#' @title Default parameter list for \code{InterpretMSSpectrum}.
#' @docType data
#' @format A data frame with 22 observations on the following 3 variables:
#'  \describe{
#'  \item{\code{ionization}}{ESI or APCI -- will influence expected peak width and precision as well as adducts.}
#'  \item{\code{ionmode}}{positive or negative -- will influence expected adducts.}
#'  \item{\code{allowed_elements}}{Passed to Rdisop in formula generation.}
#'  \item{\code{maxElements}}{Passed to Rdisop in formula generation.}
#'  \item{\code{minElements}}{Passed to Rdisop in formula generation.}
#'  \item{\code{substitutions}}{Will be deprecated in the future.}
#'  \item{\code{quick_isos}}{TRUE = via Rdisop, FALSE = via enviPat (often more correct)}
#'  \item{\code{score_cutoff}}{Specifies initial filtering step threshold per fragment. Sum Formulas with score_i < score_cutoff*max(score) will be removed.}
#'  \item{\code{neutral_loss_cutoff}}{Specifies the allowed deviation in mDa for neutral losses to be accepted from the provided neutral loss list.}
#'  }
#' @details Default parameter list used by \code{\link{InterpretMSSpectrum}}, serving also as a template 
#'     for custom lists. Basically every option which needs to be modified rarely went in here. Specific 
#'     parameter set modifications (i.e. for 'APCIpos') are provided and can be called using the character 
#'     string as a shortcut. Alternatively, a named list can be provided where all contained parameters 
#'     will receive the new specified values.
#' @usage data(param.default)
#' @keywords datasets
"param.default"