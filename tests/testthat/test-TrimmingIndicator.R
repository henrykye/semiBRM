test_that("When a matrix is given, it should return the product of trimming indicators of its columns.", {

    # data generating process
    N <- 500L
    x1 <- stats::rnorm(N)
    x2 <- stats::rt(N, df = 2)^2
    x3 <- stats::runif(N)
    x <- cbind(x1, x2, x3)

    # trimming indicator from scratch for a matrix
    alpha <- 0.025
    trimming_indicator <- rep(1L, N)
    for (k in 1:ncol(x)){

        w <- x[,k]
        cutoffs <- stats::quantile(w, c(alpha, 1-alpha))
        trimming_indicator <- trimming_indicator * ifelse(w > cutoffs[1L] & w < cutoffs[2L], 1L, 0L)
    }

    expect_equal(TrimmingIndicator(x, alpha), trimming_indicator)
})
