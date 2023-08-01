#   Revised scripts for linear, quadratic and power regression model types

#' @importFrom ggplot2 ggplot aes geom_point xlab ylab labs
base_resid_plot <- function(x, y, tit) {
  d <- as.data.frame(cbind(x, y))
  ggplot2::ggplot(d, ggplot2::aes(x, y)) +
    ggplot2::geom_point() +
    ggplot2::xlab("Fitted") +
    ggplot2::ylab("Residuals") +
    ggplot2::labs(title = tit)
}
#------------------------------------------------------------------------------
#' Plots residual plots for different regressions
#' @param d A tibble containing x (concentration) and y (response)
#' @param model_type (l)inear, (q)uadratic or (p)ower regression
#' @usage residplot(d, model_type = NULL)
#' @examples
#' data(chloromethane)
#' residplot(chloromethane)
#' residplot(chloromethane, "q")
#' @export
#' @importFrom stats lm fitted resid nls
residplot <- function(d , model_type = NULL) {
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
        base_resid_plot(stats::fitted(model),
                        stats::resid(model),
                        c("Linear"))
      },
      "q" = {
        model = stats::lm(y ~ x + I(x ^ 2))
        base_resid_plot(stats::fitted(model),
                        stats::resid(model),
                        c("Quadratic"))
      },
      "p" = {
        par <- ss.pwr(x, y)
        model <-
          minpack.lm::nlsLM(y ~ C + A * x ^ b,
                            start = list(C = par[1], A = par[2], b = par[3]))
        #        model <-
        #         stats::nls(y ~ C + A * x ^ b, # data=d,
        #                     start = list(C = par[1], A = par[2], b = par[3]))
        base_resid_plot(stats::fitted(model),
                        stats::resid(model),
                        c("Power"))
      },
      message("Check model type: ('l')inear/('q')uadratic/('p')ower")
    )
  }
}
