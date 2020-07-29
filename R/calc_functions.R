adjustcolnames <- function(d) {
    if (colnames(d[, 1]) != "x") {
        names(d)[1] <- "x"
        names(d)[2] <- "y"
    }
    return(d)
}

#' Used internally for calculating the least square regression coefficients
#' of model - assumes a linear response
#'
least_sq_reg <- function(x, y) {
    n <- length(x)
    m <- (sum(y * x) - (1/n) * sum(x) * sum(y))/(sum(x^2) - (1/n) * sum(x)^2)
    c <- 1/n * sum(y) - (m/n) * sum(x)
    return(c(slope = m, intercept = c))
}
#' Used internally to calculate residual standard deviation
#'
s_y <- function(x, y) {
    n <- length(x)
    model <- least_sq_reg(x, y)
    return(sqrt(sum((y - (model[1] * x + model[2]))^2)/(n - 2)))
    # model[1]*x + model[2] = fitted y
}
#
#' Calculates the detection limit according to Vogelgesang & Hadrich
#' @param x A vector
#' @param y A vector
#' @examples
#' dl_vogelhad(x,y)
#' @export xc
#' @export xd
dl_vogelhad <- function(x, y) {
    n <- length(x)
    ls <- least_sq_reg(x, y)
    y_crit <- ls[2] + s_y(x, y) * qt(0.95, n - 2) * sqrt(1 + 1/n + mean(x)^2/sum((x - mean(x))^2))
    x_crit <- (y_crit - ls[2])/ls[1]
    x_id <- 2 * x_crit
    return(c(xc = x_crit, xd = x_id))
}
#
#' Calculates the detection limit according to Miller & Miller
#' @param x A vector
#' @param y A vector
#' @examples
#' dl_miller(x,y)
#' @export dl
#' @export blank_dl
dl_miller <- function(x, y) {
    n <- length(x)
    ls <- least_sq_reg(x, y)
    dl <- (3 * s_y(x, y))/ls[1]
    blank_dl <- (3 * s_y(x, y) + ls[2])/ls[1]
    return(c(dl = dl, b_dl = blank_dl))
}
#' Calculates the detection limit according to Hubert & Vos using iterative calculation
#' @param x A vector
#' @param y A vector
#' @examples
#' dl_hubertvos(x,y)
#' @export dl
dl_hubertvos <- function(x, y, alpha = NULL, beta = NULL) {
    alpha <- ifelse(is.null(alpha), 0.01, alpha)
    beta <- ifelse(is.null(beta), 0.05, beta)
    n <- length(x)
    dl <- dl_vogelhad(x, y)[1]
    repeat {
        # Update dl
        new.dl <- s_y(x, y)/least_sq_reg(x, y)[1] * (qt(1 - alpha, n - 2) * sqrt(1 + 1/n +
            mean(x)^2/sum((x - mean(x))^2)) + qt(1 - beta, n - 2) * sqrt(1 + 1/n + (dl -
            mean(x))^2/sum((x - mean(x))^2)))
        conv <- sum((new.dl - dl)^2/dl^2)
        # Exit test with return() statement
        if (conv < 1e-10)
            return(dl)
        # Save interation result
        dl <- new.dl
    }
}

#' Reports detection limits ()
#' @param dat A vector containing x and y
#' @examples
#' reportDL(dat)
reportDL <- function(dat) {
    x <- dat[,1]
    y <- dat[,2]
    print(x)
    cat("Miller DL = ", round(dl_miller(x, y)[1], 1), "\n")
#    cat("Voselgang Hadrich DL = ", round(dl_vogelhad(dat[,1], dat[,2])[1], 1), "\n")
#    cat("Hubert Vos DL = ", round(dl_hubertvos(dat[,1], dat[,2])[1], 1), "\n")
}
#----------------------------------------------------
