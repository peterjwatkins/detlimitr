#------------------------ Linear detection limits ------------------------------
#
least_sq_est <- function(x, y) {
  # Used internally for calculating the least square regression coefficients
  # of model - assumes a linear response
  d.lm <- lm(y~x)
  return(c(m = coef(d.lm)[2], c = coef(d.lm)[1]))
}
s_y <- function(x, y) {
  # Used internally to calculate residual standard deviation
  #
  n <- length(x)
  ls <- least_sq_est(x, y)
  return(sqrt(sum((y - (
    ls[1] * x + ls[2]   # ls[1]*x + ls[2] = fitted y
  )) ^ 2) / (n - 2)))
}
#------------------------ Calculation ------------------------------------------
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
  cat("Linear\n")
  cat("------------------------\n")
  cat("Miller", round(dl_miller(x, y), dp), "\n")
  cat("Vogelsang-Hadrich", round(dl_vogelhad(x, y), dp), "\n")
  cat("chemCal", round(chemCal::lod((lm(
    y ~ x
  )))$x, dp), "\n")
  cat("Hubaux-Vos", round(dl_hubertvos(x, y), dp), "\n")
}
#---------------------------- Plot functions -----------------------------------
#   Original linear plot
#
##' Plots segment of calibration curve showing the estimated linear DLs.
#' @param d A tibble containing x (concentration) and y (response)
#' @usage plotlinDL(d)
#' @examples
#' data(mtbe)
#' plotlinDL(mtbe)
#' p <- plotlinDL(mtbe)
#' p + ## add additional labels
#'   ggplot2::xlab(expression(paste("Concentration (ng", " ", "g" ^ {-1}, ")")))
#' @importFrom dplyr filter
#' @importFrom chemCal lod
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_segment xlab ylab annotate
#' @export
plotlinDL <- function(d) {
  d <- adjustdf(d)
  x <- d$x
  y <- d$y
  l <- chemCal::lod(lm(y ~ x))
  model <- least_sq_est(x, y)

  xm <- dl_miller(x, y)[1]
  xv <- dl_vogelhad(x, y)[1]
  xh <- dl_hubertvos(x, y)[1]
  xc <- l$x
  ym <- model[1] * xm + model[2]
  yv <- model[1] * xv + model[2]
  yh <- model[1] * xh + model[2]
  yc <- model[1] * xc + model[2]
  xmax <- max(xm, xv, xh, xc)
  ymax <- model[1] * xmax + model[2]
  d_mod <- dplyr::filter(d, x <= xmax)
  p <-
    ggplot2::ggplot(d_mod, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = model[2],
                         slope = model[1],
                         linetype = "dotted") +
    ggplot2::xlim(c(0, 1.2 * xmax)) +
    ggplot2::ylim(c(0, 1.2 * ymax)) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        y = yv,
        xend = xv,
        yend = yv
      ),
      colour = "darkgreen",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = xv,
        y = yv,
        xend = xv,
        yend = 0
      ),
      colour = "darkgreen",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        y = ym,
        xend = xm,
        yend = ym
      ),
      colour = "blue",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = xm,
        y = ym,
        xend = xm,
        yend = 0
      ),
      colour = "blue",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        y = yh,
        xend = xh,
        yend = yh
      ),
      colour = "brown",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = xh,
        y = yh,
        xend = xh,
        yend = 0
      ),
      colour = "brown",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        y = yc,
        xend = xc,
        yend = yc
      ),
      colour = "black",
      linetype = "dashed"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = xc,
        y = yc,
        xend = xc,
        yend = 0
      ),
      colour = "black",
      linetype = "dashed"
    ) +
    ggplot2::xlab("Concentration") +
    ggplot2::ylab("Response (AU)") +
    ggplot2::annotate(
      "text",
      x = xv,
      y = 0,
      parse = TRUE,
      label = as.character(expression(italic("x")[Vogelsang - Hadrich]))
    ) +
    ggplot2::annotate(
      "text",
      x = xm,
      y = 0,
      parse = TRUE,
      label = as.character(expression(italic("x")[Miller]))
    ) +
    ggplot2::annotate(
      "text",
      x = xh,
      y = 0,
      parse = TRUE,
      label = as.character(expression(italic("x")[Hubaux - Vos]))
    ) +
    ggplot2::annotate(
      "text",
      x = xc,
      y = 0,
      parse = TRUE,
      label = as.character(expression(italic("x")[chemCal]))
    )
  return(p)
}
#-------------------------------------------------------------------------------

#' @importFrom stats lm predict coefficients
plotlinearDL <- function(x, y) {
  d.lm <- stats::lm(y ~ x)
  pred_model <-
    stats::predict(d.lm, newdata = as.data.frame(x), interval = "prediction")

  y_fit <- pred_model[, 1]
  y_lwr <- pred_model[, 2]
  y_upr <- pred_model[, 3]

  dt <-
    data.frame(cbind(x,
                     y,
                     y_fit,
                     y_lwr,
                     y_upr))
  x_dl <- dl_hubertvos(x, y)[1]
  y_dl <-
    stats::coefficients(d.lm)[1] + stats::coefficients(d.lm)[2] * x_dl

  p <- baseDLplot(dt, x_dl, y_dl, c("Linear"))

  return(p)
}
