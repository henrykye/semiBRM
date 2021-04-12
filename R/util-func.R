#' @useDynLib semiBRM, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL
#> NULL


#' Trimming indicator
#'
#' This returns an integer vector whose i-th element takes 1L if the i-th element
#' of \code{x} lies between f*100 and (1-f)*100 -th percentiles, i.e.,
#' an interior point, and 0L otherwise.
#'
#' @param x a numeric vector or matrix.
#' @param f a numeric number that indicates the cut-off percentile level.
#'
#' @return an integer vector consisting only of ones and zeros.
#'
#' If \code{x} is a matrix, interior points are searched on each column, and
#' the element-wise product of the trimming indicators of the columns are returned.
#'
#' @export
#'
#' @examples
#' n <- 500L
#' x <- rnorm(n)
#' f <- 0.025
#'
#' trimming <- TrimmingIndicator(x, f)
TrimmingIndicator <- function(x, f)
{
  if ( is.vector(x) ){

    L <- stats::quantile(x, f)
    H <- stats::quantile(x, 1-f)

    tm <- integer(length(x))
    tm[x > L & x < H] = 1L

    return(tm)
  }

  tm <- rep(1L, nrow(x))
  for ( j in seq_len(ncol(x)) ){
    tm <- tm*TrimmingIndicator(x[,j], f)
  }

  return(tm)
}


#' Nadaraya-Watson estimator with Gaussian kernel
#'
#' This estimates nonparametric conditional expectation \code{E[y|x=args]} using
#' the Nadaraya-Watson estimator with Gaussian kernel.
#'
#' @param x a numeric vector or matrix of explanatory variable(s).
#' @param y a numeric vector of outcome (or dependent) variable.
#' @param h a numeric vector indicating bandwidth size(s) for each explanatory variable in \code{x}.
#' @param args a numeric vector or matrix of arguments at which nonparametric conditional
#' expectations are evaluated.
#'
#' @details This is the Nadaraya-Watson estimator taking as arguments data \code{x} and \code{y},
#' bandwidth \code{h}, and target points \code{args}. If \code{args} are not given, it will return
#' the leave-one-out version. The dimension of \code{h} should be the same as that of \code{x}, i.e.,
#' bandwidths need to be separately specified for each explanatory variable. A popular choice of
#' bandwidth is the Silverman's rule of thumb, \code{h = sd(x)*N^(-1/r)} with \code{r = 4 + q}, where
#' \code{N} is the sample size and \code{q} is the dimension of \code{x}.
#'
#' @return a numeric vector of nonparametric estimates of conditional expectation \code{E[y|x=args]}.
#' @export
#' @examples
#' # data generating process
#' N <- 500L
#' x1 <- stats::rnorm(N)
#' x2 <- stats::rnorm(N)
#' x3 <- stats::rnorm(N)
#' x <- cbind(x1, x2, x3)
#' v <- (x1 + x2 + x3)/sqrt(3)
#' y <- v + rnorm(N)
#'
#' # univariate case
#' h_u <- stats::sd(v)*N^(-1/5)
#' vargs <- stats::quantile(v)
#'
#' yhat_univ_1 <- GaussianNadarayaWatsonEstimator(v, y, h_u, vargs)
#' yhat_univ_2 <- GaussianNadarayaWatsonEstimator(v, y, h_u) # leave-one-out version
#'
#' # multivariate case
#' h_m <- apply(matrix(stats::rnorm(150), 50, 3), 2, stats::sd)*N^(-1/7)
#' xargs <- apply(matrix(stats::rnorm(150), 50, 3), 2, stats::quantile)
#'
#' yhat_multi_1 <- GaussianNadarayaWatsonEstimator(x, y, h_m, xargs)
#' yhat_multi_2 <- GaussianNadarayaWatsonEstimator(x, y, h_m) # leave-one-out version
GaussianNadarayaWatsonEstimator <- function(x, y, h, args = NULL)
{
  stopifnot(is.vector(y))
  stopifnot(is.vector(x)|is.matrix(x))

  num.threads <- get_num_threads()

  if (num.threads == 1L){

    if (is.matrix(x)){

      if (is.null(args)){
        return( .GaussianMultivarLeaveOneOut(y, x, h) )
      }else{
        return( .GaussianMultivar(y, x, args, h) )
      }

    }else{

      if (is.null(args)){
        return( .GaussianUnivarLeaveOneOut(y, x, h) )
      }else{
        return( .GaussianUnivar(y, x, args, h) )
      }
    }

  }else if (num.threads > 1L){

    if (is.matrix(x)){

      if (is.null(args)){
        return( .GaussianMultivarLeaveOneOutOMP(y, x, h, num.threads) )
      }else{
        return( .GaussianMultivarOMP(y, x, args, h, num.threads) )
      }

    }else{

      if (is.null(args)){
        return( .GaussianUnivarLeaveOneOutOMP(y, x, h, num.threads) )
      }else{
        return( .GaussianUnivarOMP(y, x, args, h, num.threads) )
      }
    }

  }else{

    stop("The number of threads should be positive.")
  }
}


#' Marginal effects
#'
#' This computes marginal effects for a small change in a chosen continuous explanatory
#' variable as difference of semiparametric conditional probabilities.
#'
#' @param fit a fitted 'semiBRM' object.
#' @param variable an integer indicating the position of the explanatory variable of interest or
#' character of its name. For example, to compute marginal effects of the first
#' explanatory variable, \code{variable} needs to be either 1L or its name in character.
#'
#' @param delta the size of perturbation over which difference of semiparametric conditional
#' probabilities are computed. Popular choice of \code{delta} would be 1 or standard deviation of
#' the variable of interest.
#'
#' @param p.cutoffs a numeric vector of probabilities (e.g. \code{p.cutoffs = c(1/4, 2/4, 3/4)}).
#'
#' @param h a numeric of bandwidth size in the Nadaraya-Watson estimator. If not given, it will use
#' the Silverman's rule of thumb bandwidth, \code{h = sd(x)*1.06*N^(-1/5)}.
#'
#' @param trimming a logical indicating whether to trim semiparametric conditional probabilities near
#' boundaries into zeros. If \code{trimming = TRUE}, then the trimming indicator is extracted from
#' the given 'semiBRM' object.
#'
#' @details This function is designed to analyze marginal effects of a chosen explanatory variable
#' over its change by \code{delta}, where marginal effects are defined as difference between the two
#' semiparametric conditional probabilities evaluated with and without perturbation \code{delta} to
#' the variable of interest. Notice that this is designed for a continuous variable in \code{x}.
#'
#' If \code{p.cutoffs = NULL}, then it will return average marginal effects over data points. If
#' a numeric vector of probabilities is given for \code{p.cutoffs}, it will calculate quantile
#' marginal effects accordingly to the grouping by quantile values at \code{p.cutoffs}. For example,
#' if \code{p.cutoffs = c(1/4, 2/4, 3/4)}, it will return quartile marginal effects over the four
#' quartile groups.
#'
#' @return a vector or matrix of marginal effects and standard errors.
#' @export
#' @examples
#' # data generating process
#' N <- 1000L
#' X1 <- rnorm(N)
#' X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
#' X3 <- rnorm(N)^2/sqrt(2)
#' X <- cbind(X1, X2, X3)
#' beta <- c(2, 2, -1, -1)
#' V <- as.vector(cbind(X, 1)%*%beta)
#' Y <- ifelse(V >= rnorm(N), 1L, 0L)
#'
#' # parameter estimation
#' qmle <- semiBRM(x = X, y = Y, control = list(iterlim=50))
#'
#' # average marginal effects of X1
#' av_me <- MarginalEffects(qmle, variable = "X1", delta = sd(X1))
#'
#' # quantile marginal effects of X1
#' q_me <- MarginalEffects(qmle, variable = "X1", delta = sd(X1), p.cutoffs = c(1/3, 2/3))
MarginalEffects <- function(fit, variable, delta, p.cutoffs = NULL, h = NULL, trimming = FALSE)
{
  stopifnot(class(fit)=="semiBRM")
  stopifnot(is.numeric(variable)|is.integer(variable)|is.character(variable))

  if (is.numeric(variable)){
    variable <- as.integer(round(variable))
  }

  x1 <- x0 <- as.matrix(fit$model[, -1L])
  x1[, variable] <- x1[, variable] + delta

  v0 <- x0[,1L] + as.vector(x0[,-1L]%*%fit$estimate)
  v1 <- x1[,1L] + as.vector(x1[,-1L]%*%fit$estimate)

  N <- nrow(x1)
  if (is.null(h)) {h <- stats::sd(v0)*1.06*length(v0)^(-1/5)}

  prob0 <- GaussianNadarayaWatsonEstimator(v0, fit$model[,1L], h)
  prob1 <- GaussianNadarayaWatsonEstimator(v0, fit$model[,1L], h, v1)

  me <- prob1 - prob0

  if (trimming == TRUE){
    me <- me*TrimmingIndicator(v0, fit$r)
  }

  if (is.null(p.cutoffs)){

    marg.eff <- mean(me)
    std.err <- stats::sd(me)/sqrt(N)
    output <- c("Marg.Eff" = marg.eff, "Std.Err" = std.err)

  }else{

    var.cutoffs <- stats::quantile(x0[,variable], p.cutoffs)

    n_group <- length(p.cutoffs) + 1L
    group_id <- matrix(0L, nrow = N, ncol = n_group)
    group_id[, 1L] <- ifelse(x0[,variable] <= var.cutoffs[1L], 1L, 0L)
    group_id[, n_group] <- ifelse(x0[,variable] > var.cutoffs[(n_group-1L)], 1L, 0L)

    if ( n_group > 2L ){
      for (i in 1L:(length(p.cutoffs)-1L)){
        group_id[,(1L+i)] <- ifelse(var.cutoffs[i] < x0[,variable] & x0[,variable] <= var.cutoffs[(i+1)], 1L, 0L)
      }
    }

    output <- matrix(NaN, nrow = n_group, ncol = 2L)
    rownames(output) <- paste0("G", 1L:n_group)
    colnames(output) <- c("Marg.Eff", "Std.Err")

    for (g in seq_len(n_group)){

      g_N <- sum(group_id[,g])
      g_target <- group_id[,g]

      marg.eff <- sum(me[g_target==1L]) / g_N
      var <- sum( ((me-marg.eff)^2)[g_target==1L] ) / (g_N-1)
      std.err <- sqrt(var/g_N)

      output[g,] <- c(marg.eff, std.err)
    }

  }

  output
}


#' Bandwidth grid search
#'
#' This performs grid search to find the best bandwidth for the fitted semiparmetric binary
#' response model.
#'
#' @param fit a fitted 'semiBRM' object.
#' @param criterion a character indicating the criterion to be used for evaluation. It can be
#' one in c("McF.R2", "accuracy", "precision") (see Details).
#' @param n.bandwidth an integer of the number of bandwidths to be searched.
#' @param bandwidth a vector of bandwidths to be evaluated.
#'
#' @return a list with \code{bandwidth.best}: the best bandwidth found, \code{criterion.value}:
#' its criterion value, \code{criterion}: the choice of criterion, and \code{bandwith}:
#' the scoreboard for candidates.
#'
#' @details Evaluation is based on in-sample semiparametric conditional probability Pr(Y=1|X), which
#' is obtained from the leave-one-out version of the Gaussian Nadaraya-Watson estimator.
#'
#' This function first randomly generates \code{n.bandwidth} candidates from a uniform distribution
#' with range \code{[h0/3, 3*h0]}, where \code{h0} is the Silverman's rule of thumb bandwidth,
#' \code{h0 = sd(v)*1.06*N^(-1/5)}. Then, it computes 'performance' according to the given
#' \code{criterion} over the candidates and picks up the 'best' one.
#'
#' \code{criterion = "McF.R2"} employs the McFadden’s pseudo R-squared for assessment.
#'
#' For \code{criterion = "accuracy"} or \code{criterion = "precision"}, the function first predicts
#' binary response using the semiparametric conditional probability, assigning 1 if
#' the probability is greater than or equal to 0.5 and 0 otherwise. Then, it calculates the
#' confusion matrix and predictive accuracy or precision. In turn, the maximizer of the predictive
#' accuracy or precision becomes the best bandwidth.
#'
#' @export
#' @examples
#' # data generating process
#' N <- 500L
#' X1 <- rnorm(N)
#' X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
#' X3 <- rnorm(N)^2/sqrt(2)
#' X <- cbind(X1, X2, X3)
#' beta <- c(2, 2, -1, -1)
#' V <- as.vector(cbind(X, 1)%*%beta)
#' Y <- ifelse(V >= rnorm(N), 1L, 0L)
#'
#' # parameter estimation
#' qmle <- semiBRM(x = X, y = Y, control = list(iterlim=50))
#'
#' # bandwidth search
#' h1 <- BandwidthGridSearch(qmle)
#' h2 <- BandwidthGridSearch(qmle, criterion = "accuracy", bandwidth = h1$bandwidth[,1L])
#' h3 <- BandwidthGridSearch(qmle, criterion = "precision", bandwidth = h1$bandwidth[,1L])
BandwidthGridSearch <- function(fit, criterion = "McF.R2", n.bandwidth = 50, bandwidth = NULL)
{
  stopifnot(class(fit)=="semiBRM")
  stopifnot(criterion=="McF.R2"|criterion=="accuracy"|criterion=="precision")

  ## data
  y <- fit$model[,1L]
  x <- as.matrix(fit$model[,-1L])
  v <- x[,1L] + as.vector(x[,-1L]%*%fit$estimate)

  trimming_indicator <- TrimmingIndicator(v, fit$trimming.level)
  n <- length(y)
  k <- length(fit$estimate)

  ## Silverman Rule-of-Thumb and the range for grid search, suggested by Bruce E. Hansen
  h0 <- stats::sd(v)*1.06*n^(-1/5)
  h1 <- h0/3
  h2 <- 3*h0

  ## scoreboard
  if (is.null(bandwidth)){
    h_grid <- sort(stats::runif(n.bandwidth, min = h1, max = h2))
  }else{
    h_grid <- bandwidth
  }

  scoreboard <- cbind(h_grid, NaN)
  colnames(scoreboard) <- c("bandwidth", "value")
  rownames(scoreboard) <- seq_len(n.bandwidth)

  if( criterion == "McF.R2"){

    for (i in seq_along(h_grid)){

      h <- h_grid[i]
      prob <- GaussianNadarayaWatsonEstimator(v, y, h)

      prob_ <- prob[trimming_indicator == 1L]
      y_ <- y[trimming_indicator == 1L]

      logLik1 <- sum(log( c(prob_[y_ == 1L], 1-prob_[y_ == 0L]) ))
      logLik0 <- sum(y_)*log(mean(y_)) + (n-sum(y_))*log(1-mean(y_))

      mcfr2 <- 1 - logLik1/logLik0
      scoreboard[i,2L] <- mcfr2
    }

  }else if ( criterion == "accuracy" ){

    for (i in seq_along(h_grid)){

      h <- h_grid[i]
      prob <- GaussianNadarayaWatsonEstimator(v, y, h)

      prob_ <- prob[trimming_indicator == 1L]
      y_ <- y[trimming_indicator == 1L]

      y_hat <- ifelse(prob_ >= 0.5, 1L, 0L)
      conf_mat <- base::table(y_, y_hat)

      accuracy <- (conf_mat[1L,1L] + conf_mat[2L,2L]) / sum(conf_mat)
      scoreboard[i,2L] <- accuracy
    }

  }else {

    for (i in seq_along(h_grid)){

      h <- h_grid[i]
      prob <- GaussianNadarayaWatsonEstimator(v, y, h)

      prob_ <- prob[trimming_indicator == 1L]
      y_ <- y[trimming_indicator == 1L]

      y_hat <- ifelse(prob_ >= 0.5, 1L, 0L)
      conf_mat <- base::table(y_, y_hat)

      precision <- conf_mat[2L,2L] / sum(conf_mat[,2L])
      scoreboard[i,2L] <- precision
    }
  }

  h_best <- scoreboard[order(scoreboard[,2L], decreasing = TRUE),][1L,]

  list(bandwidth.best = h_best[1L],
       criterion.value = h_best[2L],
       criterion = criterion,
       bandwidth = scoreboard)
}
