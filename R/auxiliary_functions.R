#-------------------------------------------------------------------
# Auxiliary supporting functions
#
adjustdf <- function(d) {
  # Used internally for formatting the dataframe
  # Issue if the column names are not 'x' & 'y'
  if (colnames(d)[1] != "x") {
    names(d)[1] <- "x"
  }
  ## Explicitly check for y as well
  if (colnames(d)[2] != "y") {
    names(d)[2] <- "y"
  }
  return(as.data.frame(d))
}
#' @importFrom stats lm qf residuals
test_model <- function(x, y) {
  lm <- stats::lm(y ~ x)
  qm <- stats::lm(y ~ x + I(x ^ 2))
  n <- length(x)
  test_stat <-
    (sqrt(sum(stats::residuals(lm) ^ 2)) - sqrt(sum(stats::residuals(qm) ^ 2))) /
    sqrt(sum(stats::residuals(qm) ^ 2))
  f_crit <- stats::qf(0.05, 1, n - 3)

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
