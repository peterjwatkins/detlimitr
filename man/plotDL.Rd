% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_function.R
\name{plotDL}
\alias{plotDL}
\title{Plots the calibration curve.}
\usage{
plotDL(d, model_type = NULL)
}
\arguments{
\item{d}{A tibble containing x (concentration) and y (response)}

\item{model_type}{(l)inear, (q)uadratic or (p)ower regression}
}
\description{
Plots the calibration curve.
}
\examples{
data(mtbe)
plotDL(mtbe)
p <- plotDL(mtbe)
p + ## add additional labels
  ggplot2::xlab(expression(paste("Concentration (ng", " ", "g" ^ {-1}, ")")))

data(chloromethane)
plotDL(chloromethane, "q")
}
