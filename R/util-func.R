#' Interior Points Indicator
#'
#' This returns an integer vector whose i-th element takes 1L if the i-th element
#' of \code{X} is an interior point lying between \code{f} and 1-\code{f} percentiles
#' and 0L otherwise.
#'
#' @param X a vector or matrix of numeric values
#' @param f a numeric number indicating cut-offs at \code{f} and 1-\code{f} percentiles
#'
#' @return An integer vector consisting only of ones and zeros.
#'
#' If \code{X} is a matrix, interior points are searched on each column, and
#' the element-wise product of the interior point indicators among the columns are returned.
#' @export
#' @examples
#' n <- 1000L
#' X <- rnorm(n)
#' f <- 0.025
#'
#' # construct interior-point indicator from scratch
#' interior0 <- integer(n)
#' interior0[X > quantile(X, f) & X < quantile(X, 1-f)] <- 1L
#'
#' # using the function
#' interior1 <- getInteriorPoints(X, f)
#'
#' # compare outcome
#' all.equal(interior0, interior1) # TRUE
getInteriorPoints<- function(X, f)
{
  if ( is.vector(X) ){

    L <- stats::quantile(X, f)
    H <- stats::quantile(X, 1-f)

    tm <- integer(length(X))
    tm[X > L & X < H] = 1L

    return(tm)
  }

  tm <- rep(1L, nrow(X))
  for ( j in seq_len(ncol(X)) ){

    tm <- tm*getInteriorPoints(X[,j], f)
  }

  return(tm)
}
