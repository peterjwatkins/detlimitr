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

#
#' Calculates the detection limit according to Vogelgesang & Hadrich
#' @param x A vector
#' @param y A vector
#' @examples
#' dl_vogelhad(x,y)
#' @return VHdl
#' @export
dl_vogelhad <- function(x, y) {
    # Vogelgesang-Hadrich
    n <- length(x)
    ls <- least_sq_est(x, y)
    y_crit <-
        ls[2] + s_y(x, y) * qt(0.95, n - 2) *
        sqrt(1 + 1 / n + mean(x) ^ 2 / sum((x - mean(x)) ^ 2))
    x_crit <- (y_crit - ls[2]) / ls[1]
    x_id <- 2 * x_crit
    VHdl <- as.numeric(c(xc = x_crit, xd = x_id))
    return(VHdl)
}

#
#' Calculates the detection limit according to Miller & Miller
#' @param x A vector
#' @param y A vector
#' @examples
#' dl_miller(x,y)
#' @return MMdl
#' @export
dl_miller <- function(x, y) {
    #Miller - Miller
    n <- length(x)
    ls <- least_sq_est(x, y)

    dl_calc <- (3 * s_y(x, y)) / ls[1]
    dl_blank <- (3 * s_y(x, y) + ls[2]) / ls[1]
    MMdl <- as.numeric(c(dl_c = dl_calc,
              dl_b = dl_blank))
    return(MMdl)
}

#' Calculates the detection limit according to Hubert & Vos using iterative calculation
#' @param x A vector
#' @param y A vector
#' @examples
#' dl_hubertvos(x,y)
#' @return HVdl
#' @export
dl_hubertvos <- function(x, y, alpha = NULL, beta = NULL) {
    alpha <- ifelse(is.null(alpha), 0.01, alpha)
    beta <- ifelse(is.null(beta), 0.05, beta)

    n <- length(x)
    x.mean <- mean(x)
    HVdl <- dl_vogelhad(x, y)[1]

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

#' A wrapper function for reporting detection limits
#' @param d A vector containing x and y
#' @examples
#' reportDL(d)
#' @return
#' @export
reportDL <- function(d, dec_point = NULL) {
    dec_point <- ifelse(is.null(dec_point), 1, dec_point)
    d <- adjustcolnames(d)
    cat("Miller", round(dl_miller(d$x, d$y)[1],dec_point), "\n")
    cat("Vogelsang-Hadrich", round(dl_vogelhad(d$x, d$y)[1], dec_point), "\n")
    cat("Hubaux-Vos", round(dl_hubertvos(d$x, d$y)[1], dec_point), "\n")
}
#----------------------------------------------------
