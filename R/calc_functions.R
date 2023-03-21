#------------------------ Linear detection limits ------------------------------
dl_miller <- function(x, y) {
# Used internally to estimate the detection limit according to
# Miller & Miller
    n <- length(x)
    ls <- least_sq_est(x, y)
    dl_calc <- (3 * s_y(x, y)) / ls[1]
    MMdl <- as.numeric(dl_calc)
    return(MMdl)
}
# Extra code if needed
#    dl_blank <- (3 * s_y(x, y) + ls[2]) / ls[1]
#    MMdl <- as.numeric(c(dl_c = dl_calc,
#              dl_b = dl_blank))

#' @importFrom stats qt
dl_vogelhad <- function(x, y) {
    # Used internally to estimate the detection limit according to
    # Vogelgesang & Hädrich
    n <- length(x)
    ls <- least_sq_est(x, y)
    y_crit <-
        ls[2] + s_y(x, y) * qt(0.95, n - 2) *
        sqrt(1 + 1 / n + mean(x) ^ 2 / sum((x - mean(x)) ^ 2))
    x_crit <- (y_crit - ls[2]) / ls[1]
    VHdl <- as.numeric(x_crit)
    return(VHdl)
}
#   Extra code if needed
#   x_id <- 2 * x_crit
#   VHdl <- as.numeric(c(xc = x_crit, xd = x_id))

#' @importFrom stats qt
dl_hubertvos <- function(x, y, alpha = NULL, beta = NULL) {
# Used internally to estimate the detection limit acccording to
# Hubert & Vos detection limit using iterative calculation
    alpha <- ifelse(is.null(alpha), 0.01, alpha)
    beta <- ifelse(is.null(beta), 0.05, beta)
    n <- length(x)
    x.mean <- mean(x)
    HVdl <- dl_vogelhad(x, y)
    repeat {
        # Update dl
        new.dl <-
            s_y(x, y) / least_sq_est(x, y)[1] * (
                stats::qt(1 - alpha, n - 2) * sqrt(1 + 1 / n + x.mean ^ 2 / sum((x - x.mean) ^
                                                                             2)) +
                    stats::qt(1 - beta, n - 2) * sqrt(1 + 1 /
                                                   n + (HVdl - x.mean) ^ 2 / sum((x - x.mean) ^ 2))
            )
        # Compute relative error as a 2-norm.
        conv <- sum((new.dl - HVdl) ^ 2 / HVdl ^ 2)
        # Exit test with return() statement
        if (conv < 1e-10)
            return(as.numeric(HVdl))
        # Save interation result
        HVdl <- new.dl
    }
}

#--- Wrapper function for 'linear regression' functions
dl_linear <- function(x, y, dp) {
  dp <- ifelse(is.null(dp), 1, dp)
  cat("Miller", round(dl_miller(x, y), dp), "\n")
  cat("Vogelsang-Hadrich", round(dl_vogelhad(x, y), dp), "\n")
  cat("chemCal", round(chemCal::lod((lm(
    y ~ x
  )))$x, dp), "\n")
  cat("Hubaux-Vos", round(dl_hubertvos(x, y), dp), "\n")
}
#--------------  Quadratic regression calculation limits -----------------------
#' @importFrom stats lm predict
dl_quad <- function(x, y) {
  model <- stats::lm(y ~ x + I(x ^ 2))
  x_range <- as.data.frame(x)
  pred_model <-
    stats::predict(model, newdata = x_range, interval = "prediction")

  ###
  # Find x on the 'lo' curve which equals the intercept on 'hi'
  #
  # c_hi = c_lo + a_lo*x + b_lo*x^2
  # b_lo*x^2 + a_lo*x + (c_lo - c_hi) = 0
  b_lwr <- quad_beta(x, pred_model[, 2])
  b_upr <- quad_beta(x, pred_model[, 3])

  qdl <-
    quad_solver(b_lwr[3, 1], b_lwr[2, 1], b_lwr[1, 1] - b_upr[1, 1])

  return(qdl)
}


#-------------- Power regression calculation limits ----------------------------
#' @importFrom stats uniroot
dl_power <- function(x, y) {
  dt <- calc_power(x, y)
  co <- upr_lwr_pwr_coefs(dt)

  f <- function(x)
    co[2, 2] * x ^ co[3, 2] + (co[1, 1] - co[1, 2])
  pwrdl <- stats::uniroot(f, c(0, max(x)))$root
  return(pwrdl)
}

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
