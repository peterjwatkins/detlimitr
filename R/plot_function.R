#----------------------------------
#   Original linear plot
#
##' Plots segment of calibration curve showing the estimated limit for linear regression.
#' @param d A tibble containing x (concentration) and y (response)
#' @usage plotlinDL(d)
#' @examples
#' data(mtbe)
#' plotDL(mtbe)
#' p <- plotDL(mtbe)
#' p + ## add additional labels
#'   ggplot2::xlab(expression(paste("Concentration (ng", " ", "g" ^ {-1}, ")")))
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_segment xlab ylab annotate
#' @export
plotlinDL <- function(d) {
  d <- adjustcolnames(d)
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
#   Revised scripts for linear, quadratic and power regression model types

#' @importFrom ggplot2 ggplot aes geom_point xlab ylab labs
base_resid_plot <- function(x,y, tit) {
  ggplot2::ggplot(ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::xlab("Fitted") +
    ggplot2::ylab("Residuals") +
    ggplot2::labs(title = tit)
}

#' @importFrom stats lm fitted resid nls
resid_plot <- function(d , model_type = NULL) {
  d <- adjustdf(d)
  x <- d$x
  y <- d$y
  model_type <- ifelse(is.null(model_type), "l", model_type)

  if (!is.character(model_type))
    message("Model type is not known - ('l')inear/('q')uadratic/('p')ower")
  else {
    model_type <- set_model(model_type)
    switch(
      model_type,
      "l" = {
        model = stats::lm(y ~ x)
        base_resid_plot(stats::fitted(model), stats::resid(model), c("Linear"))
      },
      "q" = {
        model = stats::lm(y ~ x + I(x ^ 2))
        base_resid_plot(stats::fitted(model), stats::resid(model), c("Quadratic"))
      },
      "p" = {
        par <- sspwr(x, y)
        model <-
          stats::nls(y ~ C + A * x ^ b, # data=d,
              start = list(C = par[1], A = par[2], b = par[3]))
        base_resid_plot(stats::fitted(model), stats::resid(model), c("Power"))
      },
      message("Check model type: ('l')inear/('q')uadratic/('p')ower")
    )
  }
}
#------------------------------------------------
#' @importFrom ggplot2 ggplot aes geom_point geom_line theme geom_segment xlab ylab annotate
baseDLplot <- function(d, x_dl, y_dl, tit) {
  x <- d$x
  y <- d$y
  y_fit <- d$y_fit
  y_lwr <- d$y_lwr
  y_upr <- d$y_upr
  p <-  ggplot2::ggplot(ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::geom_line((ggplot2::aes(x, y = y_fit, col = "blue"))) +
    ggplot2::geom_line((ggplot2::aes(x, y = y_lwr, col = "red"))) +
    ggplot2::geom_line((ggplot2::aes(x, y = y_upr, col = "red"))) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::xlab("Concentration") +
    ggplot2::ylab("Response (AU)") +
    ggplot2::labs(title = tit) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        y = y_dl,
        xend = x_dl,
        yend = y_dl
      ),
      colour = "darkgreen",
      linetype = "dashed"
    )  +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = x_dl,
        y = y_dl,
        xend = x_dl,
        yend = 0
      ),
      colour = "darkgreen",
      linetype = "dashed"
    ) + ggplot2::annotate(
      "text",
      x = x_dl,
      y = 0,
      parse = TRUE,
      label = as.character(expression(italic("x")[HV]))
    )
  return(p)
}

#' @importFrom stats lm predict coefficients
plotlinearDL <- function(x, y) {
  d.lm <- stats::lm(y ~ x)
  pred_model <-
    stats::predict(d.lm, newdata = as.data.frame(x), interval = "prediction")

  y_fit <- pred_model[, 1]
  y_lwr <- pred_model[, 2]
  y_upr <- pred_model[, 3]

  dt <-
    data.frame(cbind(
      x,
      y,
      y_fit,
      y_lwr,
      y_upr
    ))
  x_dl <- dl_hubertvos(x, y)[1]
  y_dl <- stats::coefficients(d.lm)[1] + stats::coefficients(d.lm)[2] * x_dl

  p <- baseDLplot(dt, x_dl, y_dl, c("Linear"))

  return(p)
}

#' @importFrom stats lm predict
plotqDL <- function(x, y) {
  model <- stats::lm(y ~ x + I(x ^ 2))
  pred_model <-
    stats::predict(model, newdata = as.data.frame(x), interval = "prediction")

  dt <-
    data.frame(cbind(
      x,
      y,
      y_fit = pred_model[, 1],
      y_lwr = pred_model[, 2],
      y_upr = pred_model[, 3]
    ))

  x_dl <- dl_quad(x, y)
  y_dl <- quad_beta(x, pred_model[, 3])[1, 1]

  p <- baseDLplot(dt, x_dl, y_dl, c("Quadratic"))
  return(p)
}

plotpowerDL <- function(x, y) {
  dt <- calc_power(x, y)
  co <- upr_lwr_pwr_coefs(dt)

  x_dl <- dl_power(x, y)
  y_dl <- co[1, 2]

  p <- baseDLplot(dt, x_dl, y_dl, c("Power"))
  return(p)
}
#------------------------------------------------------------------------------
#  Wrapper function for DL plots
#
plotDL <- function(d, model_type = NULL) {
  d <- adjustdf(d)
  model_type <- ifelse(is.null(model_type), "l", model_type)

  if (!is.character(model_type))
    message("Model type is not known - ('l')inear/('q')uadratic/('p')ower")
  else {
    model_type <- set_model(model_type)
    switch(
      model_type,
      "l" = plotlinearDL(d$x, d$y),
      "q" = plotqDL(d$x, d$y),
      "p" = plotpowerDL(d$x, d$y),
      message("Check model type: ('l')inear/('q')uadratic/('p')ower")
    )
  }
}
