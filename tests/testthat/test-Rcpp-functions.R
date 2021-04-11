# R functions
## univariate case
RKSCE <- function(Y, X, args, h)
{
    Nargs <- length(args)
    output <- rep(NaN, Nargs)

    for ( i in seq_len(Nargs) ){
        Z <- (X-args[i])/h
        Kerh <- stats::dnorm(Z)/h
        output[i] <- sum(Y*Kerh)/sum(Kerh)
    }

    output
}

## univariate case, leave-one-out
RKSCE_leave_one_out <- function(Y, X, h)
{
    Nargs <- length(X)
    output <- rep(NaN, Nargs)

    for ( i in seq_len(Nargs) ){
        Z <- (X[-i]-X[i])/h
        Kerh <- stats::dnorm(Z)/h
        output[i] <- sum(Y[-i]*Kerh)/sum(Kerh)
    }

    output
}

## multivariate case
RKMCE <- function(Y, X, args, h)
{
    Nargs <- nrow(args)
    output <- rep(NaN, Nargs)

    for ( i in seq_len(Nargs) ){
        D <- sweep(X, MARGIN = 2, STATS = args[i,], FUN = "-", check.margin = FALSE)
        Z <- sweep(D, MARGIN = 2, STATS = h, FUN = "/", check.margin = FALSE )
        Kerh <- sweep(stats::dnorm(Z), MARGIN = 2, STATS = h, FUN = "/", check.margin = FALSE)
        w <- apply(Kerh, 1, prod)
        output[i] <- sum(w*Y)/sum(w)
    }

    output
}

## multivariate case, leave-one-out
RKMCE_leave_one_out<- function(Y, X, h)
{
    Nargs <- nrow(X)
    output <- rep(NaN, Nargs)

    for ( i in seq_len(Nargs) ){
        D <- sweep(X[-i,], MARGIN = 2, STATS = X[i,], FUN = "-", check.margin = FALSE)
        Z <- sweep(D, MARGIN = 2, STATS = h, FUN = "/", check.margin = FALSE )
        Kerh <- sweep(stats::dnorm(Z), MARGIN = 2, STATS = h, FUN = "/", check.margin = FALSE)
        w <- apply(Kerh, 1, prod)
        output[i] <- sum(w*Y[-i])/sum(w)
    }

    output
}


# Test Rcpp functions
test_that("Rcpp functions should return the same results as those written in R.", {

    ## data generating process
    N <- 500L
    x1 <- stats::rnorm(N)
    x2 <- stats::rnorm(N)
    x3 <- stats::rnorm(N)
    x <- cbind(x1, x2, x3)
    v <- (x1 + x2 + x3)/sqrt(3)
    y <- v + rnorm(N)

    ## univariate case
    h_u <- stats::sd(v)*N^(-1/5)
    vargs <- stats::quantile(v)

    yhat_univ_1 <- RKSCE(y, v, vargs, h_u)
    yhat_univ_2 <- RKSCE_leave_one_out(y, v, h_u)

    ## multivariate case
    h_m <- apply(matrix(rnorm(150), 50, 3), 2, stats::sd)*N^(-1/7)
    xargs <- apply(matrix(rnorm(150), 50, 3), 2, stats::quantile)

    yhat_multi_1 <- RKMCE(y, x, xargs, h_m)
    yhat_multi_2 <- RKMCE_leave_one_out(y, x, h_m)

    ## testing
    expect_equal(GaussianNadarayaWatsonEstimator(v, y, h_u, vargs), yhat_univ_1)
    expect_equal(GaussianNadarayaWatsonEstimator(v, y, h_u), yhat_univ_2)
    expect_equal(GaussianNadarayaWatsonEstimator(x, y, h_m, xargs), yhat_multi_1)
    expect_equal(GaussianNadarayaWatsonEstimator(x, y, h_m), yhat_multi_2)
})
