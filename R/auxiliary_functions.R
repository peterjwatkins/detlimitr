#-------------------------------------------------------------------
# Auxiliary supporting functions
#
adjustdf <- function(d) {
  # Used internally for formatting the dataframe
  # Issue if the column names are not 'x' & 'y'
  if (colnames(d)[1] != "x") {
    names(d)[1] <- "x"
    names(d)[2] <- "y"
  }
  return(as.data.frame(d))
}

test_model <- function(x, y) {
  d.lm <- lm(y ~ x)
  d.qm <- lm(y ~ x + I(x ^ 2))
  test_stat <-
    (sqrt(sum(residuals(d.lm) ^ 2)) - sqrt(sum(residuals(d.qm) ^ 2))) /
    sqrt(sum(residuals(d.qm) ^ 2))
  f_crit <- qf(0.05, 1, nrow(d) - 3)

  cat("Test statistic =", test_stat, " critical F-value", f_crit, "\n")

  if (test_stat >= f_crit)
    cat("Quadratic model recommended, not linear")
  else
    cat("Linear regression recommended")
}
#--------------------------------------------------------------------------------
set_model <- function(model_type) {
  mt <-
    ifelse(
      model_type == "l" |  model_type == "linear" ,
      "l",
      ifelse(
        model_type == "q" |  model_type == "quad",
        "q",
        ifelse(model_type == "p" | model_type == "pwr" , "p",
               "u") ## model type is 'u'nkwown
      )
    )
  return(mt)
}
#-------------------------------------------------------------------------------
#     Detection limits
#-------------------------------------------------------------------------------
#
#      Linear
#
least_sq_est <- function(x, y) {
  # Used internally for calculating the least square regression coefficients
  # of model - assumes a linear response
  # Included for pedagical reasons, and completeness
  n <- length(x)
  m <- (sum(y * x) - 1 / n * sum(x) * sum(y)) /
    (sum(x ^ 2) - 1 / n * sum(x) ^ 2)
  c <- 1 / n * sum(y) - m / n * sum(x)
  return(c(m = m, c = c))
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
#
#      Quadratic
#
quad_beta <- function(x, y) {
  X <- as.matrix(cbind(
    x0 = rep(1, length(x)),
    x1 = x,
    x2 = x ^ 2
  ))
  Y <- as.matrix(y)

  beta = solve(t(X) %*% X) %*% (t(X) %*% Y)
  return(beta)
}
quad_solver <- function(a, b, c) {
  x1 <- (-b + sqrt((b ^ 2 - 4 * a * c))) / (2 * a)
  x2 <- (-b - sqrt((b ^ 2 - 4 * a * c))) / (2 * a)
  if (x1 >= 0)
    x <- x1
  else
    x <- x2
  return(x)
}
#
#      Power
#
#define function to minimize residual sum of squares
min_rss <- function(data, par) {
  with(data, sum((par[1] + par[2] * x ^ par[3] - y) ^ 2))
}

sspwr <- function(x, y) {
  p.lm <- stats::lm(log1p(y) ~ log1p(x))
  df <- as.data.frame(cbind(x, y))
  res <- stats::optim(
    par = c(y[1], 10 ^ coef(p.lm)[1], 10 ^ coef(p.lm)[2]),
    fn = min_rss,
    data = df
  )
  # c(C,A,b)
  return(c(res$par[1], res$par[2], res$par[3]))
}
calc_power <- function(x, y) {
  d <- as.data.frame(cbind(x, y))

  par <- sspwr(x, y)
  model <-
    stats::nls(y ~ C + A * x ^ b, # data=d,
               start = list(C = par[1], A = par[2], b = par[3]))

  pred <- ggtrendline::predFit(
    model ,
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

upr_lwr_pwr_coefs <- function(df) {
  par <- sspwr(df$x, df$y)

  dt <- tibble::as_tibble(cbind(x = df$x, y = df$y_lwr))
  dt <- dplyr::filter(dt, x > 0)
  lwr <- stats::nls(y ~ C + A * x ^ b,
                    data = dt,
                    start = list(C = par[1], A = par[2], b = par[3]))


  dt <- tibble::as_tibble(cbind(x = df$x, y = df$y_upr))
  dt <- dplyr::filter(dt, x > 0)
  upr <- stats::nls(y ~ C + A * x ^ b,
                    data = dt,
                    start = list(C = par[1], A = par[2], b = par[3]))

  return(cbind(coef(lwr), coef(upr)))
}
#-------------------------------------------------------------------
