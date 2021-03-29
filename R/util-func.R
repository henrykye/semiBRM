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
#' n <- 1000L
#' x <- rnorm(n)
#' f <- 0.025
#'
#' # construct a trimming indicator by hand
#' ti0 <- integer(n)
#' ti0[x > quantile(x, f) & x < quantile(x, 1-f)] <- 1L
#'
#' # use the Trimming Indicator function
#' ti1 <- TrimmingIndicator(x, f)
#'
#' # compare outcome
#' all.equal(ti0, ti1) # TRUE
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
#' This estimates nonparametrically conditional expectation \code{E[y|x=args]} using
#' the Nadaraya-Watson estimator. If target points \code{args} are not given, it will return
#' the leave-one-out version.
#'
#' @param x a numeric vector or matrix of explanatory variable(s).
#' @param y a numeric vector of outcome variable.
#' @param h a numeric vector indicating bandwidth size(s) for each explanatory variable in \code{x}.
#' @param args a numeric vector or matrix of arguments at which nonparametric conditional
#' expectations are evaluated.
#'
#' @return a numeric vector of nonparametric estimates of conditional expectation \code{E[y|x=args]}.
#' @export
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
#' This computes marginal effects with respect to a chosen explanatory variable
#' as difference of semiparametric conditional probabilities.
#'
#' @param fit a fitted 'semiBRM' object.
#' @param variable an integer indicating the position of the explanatory variable of interest or
#' character of its name. For example, to compute marginal effects of the first
#' explanatory variable, \code{variable} needs to be either 1L or its name in character.
#'
#' @param delta the size of perturbation over which difference of semiparametric conditional
#' probabilities are computed. Populara choice of \code{delta} would be 1 or standard deviation of
#' the variable of interest.
#'
#' @param p.cutoffs a numeric vector of probabilities (e.g. \code{p.cutoffs = c(1/4, 2/4, 3/4)}).
#'
#' @param trimming a logical indicating whether to trim semiparametric conditional probabilities near
#' boundaries into zeros. If \code{trimming = TRUE}, then the trimming indicator is extracted from
#' the given 'semiBRM' object.
#'
#' @details This function is designed to analyze marginal effects of a chosen explanatory variable
#' over its change by \code{delta}, where marginal effects are defined as difference between the two
#' semiparametric conditional probabilites evaluated with and without perturbation \code{delta} to
#' the variable of interest.
#'
#' If \code{p.cutoffs = NULL}, then it will return average marginal effects over data points. If
#' a numeric vector of probabilities is given for \code{p.cutoffs}, it will calculate quantile
#' marginal effects accordingly to the grouping by quantile values at \code{p.cutoffs}. For example,
#' if \code{p.cutoffs = c(1/4, 2/4, 3/4)}, it will return quartile marginal effects over the four
#' quartile groups.
#'
#' @return a vector or matrix of marginal effects, standard errors, t-values, and p-values.
#' @export
MarginalEffects <- function(fit, variable, delta, p.cutoffs = NULL, trimming = FALSE)
{
  stopifnot(is.numeric(variable)|is.integer(variable)|is.character(variable))

  if (is.numeric(variable)){
    variable <- as.integer(round(variable))
  }

  x1 <- x0 <- as.matrix(fit$model[, -1L])
  x1[, variable] <- x1[, variable] + delta

  v0 <- x0[,1L] + as.vector(x0[,-1L]%*%fit$estimate)
  v1 <- x1[,1L] + as.vector(x1[,-1L]%*%fit$estimate)

  N <- nrow(x1)
  h <- stats::sd(v0)*N^(-1/5)

  prob0 <- GaussianNadarayaWatsonEstimator(v0, fit$model[,1L], h)
  prob1 <- GaussianNadarayaWatsonEstimator(v0, fit$model[,1L], h, v1)

  me <- prob1 - prob0

  if (trimming == TRUE){
    me <- me*TrimmingIndicator(v0, fit$r)
  }

  if (is.null(p.cutoffs)){

    marg.eff <- mean(me)
    std.err <- stats::sd(me)/sqrt(N)

    t.value <- round(marg.eff/std.err, 2)
    p.value <- round(2*stats::pnorm(-abs(t.value)), 4)

    output <- c("Marg.Eff" = marg.eff, "Std.Err" = std.err, "t.value" = t.value, "p.value" = p.value)

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

    output <- matrix(NaN, nrow = n_group, ncol = 4L)
    rownames(output) <- paste0("G", 1L:n_group)
    colnames(output) <- c("Marg.Eff", "Std.Err", "t.value", "p.value")

    for (g in seq_len(n_group)){

      g_N <- sum(group_id[,g])
      g_target <- group_id[,g]

      marg.eff <- sum(me[g_target==1L]) / g_N
      var <- sum( ((me-marg.eff)^2)[g_target==1L] ) / (g_N-1)
      std.err <- sqrt(var/g_N)

      t.value <- round(marg.eff/std.err, 2)
      p.value <- round(2*stats::pnorm(-abs(t.value)), 6)

      output[g,] <- c(marg.eff, std.err, t.value, p.value)
    }

  }

  output
}
