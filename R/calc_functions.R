#--------------- Main detection limit function----------------------------------
#' Summarises detection limits
#' @description A generic function summarising the detection limits (DLs) according
#' to linear, quadratic or power regression. For linear regression, the DLs are
#' estimated using four approaches;i) Miller and Miller, ii) Vogelgesang and Hädrich,
#' iii) Hubert & Vos, and iv) the R \emph{chemCal} package. The first two DLs are
#' directly estimated using internal functions, \emph{dl_miller} and \emph{dl_vogelhad}
#' while the third DL is estimated iteratively using \emph{dl_hubertvos}.
#' The default regression is assumed to be linear. However, the DLs can also be
#' estimated for either a quadratic ("q") or power ("p") regression.
#' For quadratic and power regression, the DLs are estimated using the intersection
#' of the upper prediction limit with the *x* = 0 line which then intersects
#' the lower prediction line. The associated *x* value is taken as the estimated detection limit.
#' See v) for further details on this approach. The user is advised
#' though to verify the best regression approach by assessing the residual
#' plots - see \emph{resid_plot} for further details.
#'
#' The default number of single decimal point is 1 but this can be changed by the user.
#' @param d A tibble containing x (concentration) and y (response).
#' @param model_type for regression: (l)inear, (q)uadratic or (p)ower.
#' @param dp Number of decimal points
#' @usage calcDL(d, model_type = NULL, dp = NULL)
#' @source {
#' i) J.C. Miller and J.N. Miller (1993), "Statistics for Analytical Chemistry",
#' 3rd ed., Prentice-Hall.
#' ii) J. Vogelgesang and J. Hädrich (1998), Accred. Qual. Assur., 3:242-255.
#' iii) A. Hubaux and G. Vos (1970), Anal. Chem., 42:849-855
#' & D.T. O'Neill, E.A. Rochette and P.J. Ramsay, (2002), Anal. Chem., 74:5907-5911
#' iv) J. Ranke, (2018), chemCal: Calibration Functions for Analytical Chemistry,
#' https://CRAN.R-project.org/package=chemCal
#' v) D. Coleman and L. Vanatta (2009). American Laboratory, Statistics in Analytical Chemistry: Part 34 - Detection Limit Summary.
#' }
#' @examples
#' data(mtbe)
#' calcDL(mtbe) 	  # single decimal point, with default to linear regression
#' calcDL(mtbe, 3)  # error as regression type not defined
#' calcDL(mtbe, dp = 3)
#'
#' data(chloromethane)
#' calcDL(chloromethane, dp = 2)
#' calcDL(chloromethane, "q", dp = 2) # quadratic regression
#' @export
calcDL <- function(d, model_type = NULL, dp = NULL) {
  d <- adjustdf(d)
  dp <- ifelse(is.null(dp), 1, dp)

  model_type <- ifelse(is.null(model_type), "l", model_type)

  if (is.numeric(model_type))
    message("Model type is numeric: use 'dp =' to set number of points")
  else
    if (!is.character(model_type))
      message("Model type is not known - ('l')inear/('q')uadratic/('p')ower")
  else {
    model_type <- set_model(model_type)
    switch(
      model_type,
      "l" = dl_linear(d$x, d$y, dp),
      "q" = cat("Vogelsang-Hadrich", round(dl_quad(d$x, d$y), dp), "\n"),
      "p" = cat("Vogelsang-Hadrich", round(dl_power(d$x, d$y), dp), "\n"),
      message("Check model type: ('l')inear/('q')uadratic/('p')ower")
    )
  }
}

#----------------------------------------------------
