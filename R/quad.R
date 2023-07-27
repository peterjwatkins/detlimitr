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

#------------------------ Plot -------------------------------------------------
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



