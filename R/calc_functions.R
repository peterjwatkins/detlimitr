adjustcolnames <- function(d) {
# Used internally for formatting the dataframe
# Issue if the column names are not 'x' & 'y'
    if (colnames(d)[1] != "x") {
        names(d)[1] <- "x"
        names(d)[2] <- "y"
    }
    return(d)
}
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
dl_vogelhad <- function(x, y) {
# Used internally to estimate the detection limit acccording to
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
dl_miller <- function(x, y) {
# Used internally to estimate the detection limit acccording to
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
                qt(1 - alpha, n - 2) * sqrt(1 + 1 / n + x.mean ^ 2 / sum((x - x.mean) ^
                                                                             2)) +
                    qt(1 - beta, n - 2) * sqrt(1 + 1 /
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
#' Summarises detection limits
#' @description A generic function summarising the detection limits (DLs) estimated using three approaches;
#' i) Miller and Miller, ii) Vogelgesang and Hädrich, and iii) Hubert & Vos.
#' The first two DLs are directly estimated using internal functions, \emph{dl_miller} and \emph{dl_vogelhad}
#' while the third DL is estimated iteratively using \emph{dl_hubertvos}.
#' By default, a single decimal point is shown but this can be changed by the user.
#' @param d A tibble containing x (concentration) and y (response).
#' @param dp Number of decimal points
#' @usage summaryDL(d, dp = 1)
#' @source {
#' i) J.C. Miller and J.N. Miller (1993), "Statistics for Analytical Chemistry",
#' 3rd ed., Prentice-Hall.
#' ii) J. Vogelgesang and J. Hädrich (1998), Accred. Qual. Assur., 3:242-255.
#' iii) A. Hubaux and G. Vos (1970), Anal. Chem., 42:849-855
#' & D.T. O'Neill, E.A. Rochette and P.J. Ramsay, (2002), Anal. Chem., 74:5907-5911
#' }
#' @examples
#' data(mtbe)
#' summaryDL(mtbe) 	    #single decimal point
#' summaryDL(mtbe, 3)	#three decimal points
#' @return
#' @export
summaryDL <- function(d, dp = NULL) {
    dp <- ifelse(is.null(dp), 1, dp)
    d <- adjustcolnames(d)
    cat("Miller", round(dl_miller(d$x, d$y),dp), "\n")
    cat("Vogelsang-Hadrich", round(dl_vogelhad(d$x, d$y), dp), "\n")
    cat("Hubaux-Vos", round(dl_hubertvos(d$x, d$y), dp), "\n")
}
#----------------------------------------------------
