#' Tabulate detection limits
#' @description A function which tabulates the detection limits estimated using four approaches;
#' i) Miller and Miller, ii) Vogelgesang and Hädrich, iii) Hubert & Vos, and iv) the R \emph{chemCal} package.
#' By default, a single decimal point is shown but this can be changed by the user.
#' @param d A tibble containing x (concentration) and y (response).
#' @param dp Number of decimal points
#' @usage tabulateDL(d, dp = NULL)
#' @examples
#' data(mtbe)
#' tabulateDL(mtbe) 	 #single decimal point
#' tabulateDL(mtbe, 3)	#three decimal points
#' @importFrom dplyr tibble
#' @importFrom gt gt
#' @importFrom chemCal lod
#' @importFrom stats lm
#' @export
tabulateDL <- function(d, dp = NULL) {
  dp <- ifelse(is.null(dp), 1, dp)
  d <- adjustcolnames(d)
  l <- (stats::lm(y ~ x, data = d))

  DL_tbl <-
    dplyr::tibble(
      Method = c(
        "Vogelsang - H\u00E4drich",
        "Miller - Miller",
        "chemCal",
        "Hubaux - Vos"
      ),
      DL = c(
        round(dl_vogelhad(d$x, d$y), dp),
        round(dl_miller(d$x, d$y), dp),
        round(chemCal::lod(l)$x, dp),
        round(dl_hubertvos(d$x, d$y), dp)
      )
    )
  DL_tab <- gt::gt(data = DL_tbl)
  print(DL_tab)
}
