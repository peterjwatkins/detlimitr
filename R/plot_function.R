##' Plots segment of calibration curve showing the estimated limits.
#' @param d A tibble containing x (concentration) and y (response)
#' @usage plotDL(d)
#' @examples
#' data(mtbe)
#' plotDL(mtbe)
#' p <- plotDL(mtbe)
#' p + ## add additional labels
#'   ggplot2::xlab(expression(paste("Concentration (ng", " ", "g" ^ {-1}, ")")))
#' @export
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_segment xlab ylab annotate
plotDL <- function(d) {
  d <- adjustcolnames(d)
  x <- d$x
  y <- d$y
  model <- least_sq_est(x, y)

  xm <- dl_miller(x, y)[1]
  xv <- dl_vogelhad(x, y)[1]
  xh <- dl_hubertvos(x, y)[1]
  ym <- model[1] * xm + model[2]
  yv <- model[1] * xv + model[2]
  yh <- model[1] * xh + model[2]
  xmax <- max(xm, xv, xh)
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
    ggplot2::geom_segment(ggplot2::aes(
      x = 0,
      y = yv,
      xend = xv,
      yend = yv
    ),
    colour = "darkgreen",
    linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(
      x = xv,
      y = yv,
      xend = xv,
      yend = 0
    ),
    colour = "darkgreen",
    linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(
      x = 0,
      y = ym,
      xend = xm,
      yend = ym
    ),
    colour = "blue",
    linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(
      x = xm,
      y = ym,
      xend = xm,
      yend = 0
    ),
    colour = "blue",
    linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(
      x = 0,
      y = yh,
      xend = xh,
      yend = yh
    ),
    colour = "brown",
    linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(
      x = xh,
      y = yh,
      xend = xh,
      yend = 0
    ),
    colour = "brown",
    linetype = "dashed") +
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
    )
  return(p)
}
