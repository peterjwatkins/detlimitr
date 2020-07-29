#' Plots segment of calibration curve showing the position of the different DLs
#' @param dat A vector containing x and y
#' @examples
#' show_limits(dat)
#' @export
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot aes geom_point geom_abline geom_segment xlab ylab annotate
show_limits <- function(dat) {
  x <- dat$x
  y <- dat$y
  ls <- least_sq_reg(x, y)
  # ls[1] = slope ls[2] = intercept
  x.seq <- seq(0, max(x), 0.1)
  y.fit <- ls[1] * x.seq + ls[2]
  xm <- dl_miller(x, y)[1]
  ym <- ls[1] * xm + ls[2]
  xv <- dl_vogelhad(x, y)[1]
  yv <- ls[1] * xv + ls[2]
  xh <- dl_hubertvos(x, y)[1]
  yh <- ls[1] * xh + ls[2]
  xmax <- max(xm, xv, xh)
  ymax <- ls[1] * xmax + ls[2]

  d <- dplyr::filter(dat, x <= xmax)

  p <- ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::geom_abline(intercept = ls[2], slope = ls[1], linetype = "dotted") +
    ggplot2::xlim(c(0, 1.25 * xmax)) + ggplot2::ylim(c(0,1.25 * ymax)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, y = yv, xend = xv, yend = yv),
                          colour = "darkgreen", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = xv,y = yv, xend = xv, yend = 0),
                          colour = "darkgreen", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = 0,y = ym, xend = xm, yend = ym),
                          colour = "blue", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = xm,y = ym, xend = xm, yend = 0),
                          colour = "blue", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = 0,y = yh, xend = xh, yend = yh),
                          colour = "brown", linetype = "dashed") +
    ggplot2::geom_segment(ggplot2::aes(x = xh,y = yh, xend = xh, yend = 0),
                          colour = "brown", linetype = "dashed") +
    ggplot2::xlab("Concentration") +
    ggplot2::ylab("Response (AU)") +
    ggplot2::annotate("text", x = xv, y = 0, parse = TRUE,
        label = as.character(expression(italic("x")[Vogelsang - Hadrich]))) +
    ggplot2::annotate("text", x = xm, y = 0, parse = TRUE,
        label = as.character(expression(italic("x")[Miller]))) +
    ggplot2::annotate("text", x = xh, y = 0, parse = TRUE,
        label = as.character(expression(italic("x")[Hubaux - Vos])))
  return(p)
}
