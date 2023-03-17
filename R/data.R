#' Calibration data from Table 6, Vogelsang et al (1998).
#'
#' A sample dataset of concentration (μg/kg) and measured response (area).
#'
#' @format A data frame with 12 rows and 2 variables:
#' \describe{
#'   \item{concentration}{μg/kg}
#'   \item{area}{arbitrary units}
#' }
#' @usage data(vhtab6)
#' @source {J. Vogelgesang and J. Hädrich (1998), Accred. Qual. Assur., 3:242-255.}
"vhtab6"

#' Calibration dataset for MTBE in water using GC-MS.
#'
#' A dataset containing methyl tert-butyl ether (MTBE) concentration (μg/L, x) and measured response (area, y).
#'
#' @format A data frame with 18 rows and 2 variables:
#' \describe{
#'   \item{x}{μg/L}
#'   \item{y}{arbitrary units}
#' }
#' @usage data(mtbe)
#' @source {D.T. O'Neill, E.A. Rochette and P.J. Ramsay, (2002), Anal. Chem., 74:5907-5911.}
"mtbe"

#' Paraquat calibration data using multiple square wave voltammetry.
#'
#' A dataset containing concentration (μmol/L, x) and measured response (μA X -1, y).
#'
#' @format A data frame with 9 rows and 2 variables:
#' \describe{
#'   \item{x}{μmol/L}
#'   \item{y}{μA (X -1)}
#' }
#' @usage data(paraquat)
#' @source {O.B. da Silva and S.A.S. Machado, (2012), Anal. Methods, 4:2348-2354}
"paraquat"

#' Chloromethane calibration data using GC-MS
#'
#' Chloromethane (CH3Cl) calibration data using GC-MS with fluorobenzene as an internal standard (IS)
#'
#' A dataset containing chloromethane concentration (μg/L, x) and the ratio of the measured peak
#' area for chloromethane to that of the internal standard, flurobenzene.
#'
#' @format A data frame with 90 rows and 2 variables:
#' \describe{
#'   \item{Concentration}{μg/L}
#'   \item{Response}{Arbitrary unit}
#' }
#' @usage data(chloromethane)
#' @source {I. Lavagnini and F. Magno, (2007), Mass Spec. Rev., 26:1-18}
"chloromethane"

#' HPLC calibration data for five analytes.
#'
#' Calibration data for five analytes (ibuprofen, genisten, biochanin,
#' pseudoephedrine and sodium nitrate) using different instrumental
#' conditions.
#'
#' @format A data frame with 70 rows and three variables:
#' \describe{
#'    \item{Analyte}{}
#'    \item{Concentration}{Different units}
#'    \item{Area}{Arbitrary unit}
#' }
#' @usage data(hplccal)
#' @source {L. Kirkup and M. Mulholland, (2004), J. Chrom. A, 1029:1-11}
"hplccal"
