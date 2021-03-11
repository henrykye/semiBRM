#' Quasi maximum likelihood estimation
#'
#' This implements quasi maximum likelihood estimation using "maxLik::maxLik" to find identifiable parameters.
#' This runs the estimation twice, where the first estimation is the "pilot" version and the second
#' is the main one with interior points searched based on the "index" obtained from the pilot version.
#'
#' @param X a matrix of explanatory variables.
#' @param Y an integer vector of binary response outcome
#' @param r a numeric number that controls the size of Silverman's rule-of-thumb bandwidth, h = sd(x)*N^(-r)
#' @param cutoff.level a numeric number to generate an interior-point indicator
#' @param num.threads an integer number indicating the number of threads to be used. The default is
#' parallel::detectCores()-1L
#' @param ... further arguments in maxLik::maxLik
#'
#' @return object of class 'maxLik'
#' @export
semiBRM <- function(X, Y, r = 1/7.1, cutoff.level = 0.025, num.threads = parallel::detectCores()-1L, ...)
{
  # setup
  N <- length(Y)
  # Xs <- sweep(X, 2, matrixStats::colSds(X), "/") # rescaling

  # log likelihood function
  logLike <- function(parms)
  {
    V <- X[,1L] + as.vector(X[,-1L]%*%parms)
    H <- stats::sd(V)*N^(-r)

    prob <- GaussianKerNonparLeaveOneOut(Y, V, H, num.threads)

    if ( is.null(interior.points) ){
      loglik <- log( c(prob[Y==1L], 1-prob[Y==0L]) )
      return( sum(loglik) )
    }

    prob_ <- prob[interior.points == 1L]
    Y_ <- Y[interior.points == 1L]

    loglik <- log( c(prob_[Y_ == 1L], 1-prob_[Y_ == 0L]) )
    return( sum(loglik) )
  }

  # starting values
  bols <- stats::coefficients(RcppArmadillo::fastLm(cbind(X, 1), Y))
  start <- (bols[2L:(length(bols)-1L)]/bols[1L])

  # pilot estimation
  interior.points <- semiBRM::getInteriorPoints(X, cutoff.level)
  quasi_mle_0 <- maxLik::maxLik(logLike, start = start, method = "BFGS", ...)

  # main estimation
  Vhat <- X[,1L] + as.vector(X[,-1L]%*%quasi_mle_0$estimate)
  interior.points <- semiBRM::getInteriorPoints(Vhat, cutoff.level)
  quasi_mle_1 <- maxLik::maxLik(logLike, start = quasi_mle_0$estimate, method = "BFGS", ...)

  # return
  return(quasi_mle_1)
}
