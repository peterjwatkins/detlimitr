#
#      Power
#
#' @importFrom stats lm coef
# Power function: y = C + A*x^b
ss.pwr <- function(x, y) {
  y.d <- y - min(y)
  p.lm <- stats::lm(log10(y.d+1) ~ log10(x+1))
  A <- 10^coef(p.lm)[1]
  b <- coef(p.lm)[2]
  C <- min(y) + A * x[1]^b
  # structure of returned data ->c(C,A,b)
  params <- c(C, A, b)
  names(params) <- c("","","")
  return(params)
}

#' @importFrom minpack.lm nlsLM
#' @importFrom ggtrendline predFit
calc_power <- function(x, y) {
  d <- as.data.frame(cbind(x, y))

  par <- ss.pwr(x, y)

  m1 <-
    minpack.lm::nlsLM(y ~ C + A * x ^ b,
                      start = list(C = par[1], A = par[2], b = par[3]))
#  model <-
#    stats::nls(y ~ C + A * x ^ b,
#               start = list(C = par[1], A = par[2], b = par[3]))

  pred <- ggtrendline::predFit(
    m1,
    data.frame(x = x),
    se.fit = TRUE,
    level = 0.95,
    interval = "prediction"
  )
  y_fit <- pred$fit[, 1]
  y_lwr <- pred$fit[, 2]
  y_upr <- pred$fit[, 3]

  d <- as.data.frame(cbind(x = x, y = y, y_fit, y_lwr, y_upr))
  return(d)
}
#' @importFrom stats nls coef
#' @importFrom tibble as_tibble
#' @importFrom dplyr filter


upr_lwr_pwr_coefs <- function(df) {
  xt <- df$x
  yt <- df$y
  y_lwr <- df$y_lwr
  y_upr <- df$y_upr

  par <- ss.pwr(xt, yt)

  dt <- tibble::as_tibble(cbind(xt, y_lwr))
  dt <- dplyr::filter(dt, xt > 0)

  lwr <- minpack.lm::nlsLM(y_lwr ~ C + A * xt ^ b, data = dt,
                           start = list(C = par[1], A = par[2], b = par[3]))

#  lwr <- stats::nls(y ~ C + A * x ^ b,
#                    data = dt,
#                    start = list(C = par[1], A = par[2], b = par[3]))

  dt <- tibble::as_tibble(cbind(xt, y_upr))
  dt <- dplyr::filter(dt, xt > 0)
  upr <- minpack.lm::nlsLM(y_upr ~ C + A * xt ^ b, data = dt,
                    start = list(C = par[1], A = par[2], b = par[3]))

#  upr <- stats::nls(y ~ C + A * x ^ b,
##                    data = dt,
#                    start = list(C = par[1], A = par[2], b = par[3]))

  return(cbind(stats::coef(lwr), stats::coef(upr)))
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

#-------------- Plot -----------------------------------------------------------
plotpowerDL <- function(x, y) {
  dt <- calc_power(x, y)
  co <- upr_lwr_pwr_coefs(dt)

  x_dl <- dl_power(x, y)
  y_dl <- co[1, 2]

  p <- baseDLplot(dt, x_dl, y_dl, c("Power"))
  return(p)
}







