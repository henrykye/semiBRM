#' @importFrom Rdpack reprompt
NULL
#> NULL


#' Semiparametric Binary Response Model
#'
#' Implementation of semiparametric binary response models theorized in
#' \insertCite{klein1993efficient;textual}{semiBRM}.
#'
#' @details
#' This package offers an implementation of semiparametric binary response models, theorized in
#' a seminal work in the semiparametric econometrics literature, \insertCite{klein1993efficient;textual}{semiBRM}.
#'
#' Compared to other related packages that help run non-/semi-parametric analysis,
#' this relfects the econometrician's perspective to the best, taking conditions
#' for asymptotic properties seriously. For instance, the default choice of bandwidth size in
#' the Nadaraya-Watson estimator meets the conditions for square-root N consistency of the coefficient
#' estimator.
#'
#' In turn, this package will be useful in conducing econometric analysis on binary response models.
#' For example, this package offers computation of marginal effects, which are often the most
#' important quantity of interest in binary response models in the economics context.
#'
#' This is built upon \code{\link[Rcpp:Rcpp-package]{Rcpp}} \insertCite{eddelbuettel2011rcpp}{semiBRM} along with
#' \code{OpenMP}, speeding up computation of nonparametric conditional expectation over data points.
#' In author's opinion, the only significant disadvantage of this semiparametric approach to binary
#' choice models over 'standard' ones such as Probit and Logistic, particularly in the cross-sectional
#' setup, is high computation costs. In this package, this limitation is meaningfully overcome by
#' taking advantage of multithreading via \code{OpenMP} in \code{\link[Rcpp:Rcpp-package]{Rcpp}}.
#'
#' The econometric theory underlying this package mostly comes from \insertCite{klein1993efficient;textual}{semiBRM}.
#' However, a few important parts are based on lectures of Prof. Klein at Rutgers University who
#' wrote that paper. Since \insertCite{klein1993efficient;textual}{semiBRM} was published,
#' he made several improvements in asymptotic theories. For example, the 'current version' of
#' asymptotics no longer requires the use of 'higher-order' kernel in the Nadaraya-Watson estimator,
#' which was asked in the original paper.
#'
#' @author Hyungil Kye (Henry) `<`\email{hkye@@economics.rutgers.edu}`>`
#'
#' @references
#' \insertAllCited{}
#'
#' @docType package
#' @name semiBRM-package
NULL
#> NULL


# constant variable controlling the number of the threads for parallelizing computation
.pkgglobalenv <- new.env(parent=emptyenv())
.pkgglobalenv$num.threads <- parallel::detectCores()-1L


#' Set the number of threads
#'
#' This sets the number of threads to parallelize computation of nonparametric conditional expectation.
#'
#' @param x an integer indicating the number of the threads put in place for parallelization.
#'
#' @export
set_num_threads <- function(x = parallel::detectCores()-1L)
{
    assign("num.threads", as.integer(round(x)), envir = .pkgglobalenv)
    return(invisible(NULL))
}


#' Get the number of threads
#'
#' This returns the number of threads set to parallelize computation of nonparametric conditional expectation.
#'
#' @return an integer indicating the number of the threads set for parallelization.
#' @export
get_num_threads <- function()
{
    .pkgglobalenv$num.threads
}


#' Semiparametric binary response model: Parameter estimation
#'
#' This implements quasi maximum likelihood estimation using "\code{\link[maxLik:maxLik]{maxLik::maxLik}}"
#' to estimate an identifiable  set of parameters. This runs estimation twice, where the first estimation
#' is the "pilot" version and the second is the main one with trimming indicator obtained from the "index"
#' of the pilot version in place.
#'
#' @param x a numeric matrix of explanatory variables.
#' @param y a vector of integer, numeric, or factor of binary response outcomes, taking either 1 or 0 only.
#' @param formula a formula describing the model to be fitted..
#' @param data a data.frame containing variables in \code{formula}.
#' @param r a numeric number that controls the size of Silverman's rule-of-thumb bandwidth, \code{h = sd(x)*N^(-r)}.
#' @param tau a numeric indicating cut-off levels for trimming in \code{\link[=TrimmingIndicator]{TrimmingIndicator}(X,f)},
#' which assigns 1L to the values in \code{X} lying between \code{tau}*100 and (1-\code{tau})*100 -th percentiles and 0L to
#' those outside this range.
#'
#' @param ... further arguments in "\code{\link[maxLik:maxLik]{maxLik::maxLik}}" such as \code{control}.
#'
#' @name semiBRM
#'
#' @return object of class 'semiBRM' similar to that of 'maxLik'
#' \itemize{
#'   \item \code{estimate}: estimated parameter values.
#'   \item \code{log.likelihood}: log likelihood at the estimates.
#'   \item \code{gradient}: a gradient vector at the estimates.
#'   \item \code{hessian}: a hessian matrix at the estimates.
#'   \item \code{code}: return code as detailed in \code{\link[maxLik:maxLik]{maxLik::maxLik}}.
#'   \item \code{message}: a short message describing return code.
#'   \item \code{iter}: the number of iterations performed for numerical optimization.
#'   \item \code{control}: the optimization control parameters as detailed in \code{\link[maxLik:maxLik]{maxLik::maxLik}}.
#'   \item \code{model}: the model frame.
#'   \item \code{r}: the bandwidth parameter for Silverman's rule-of-thumb bandwidth.
#'   \item \code{trimming.level}: the trimming cutoff level, which is \code{tau} in function argument.
#'   \item \code{call}: the matched call.
#'   \item \code{formula}: the formula entered for estimation.
#' }
#'
#' @details
#' This can take as arguments either matrix \code{x} and vector \code{y} or \code{formula} and \code{data} to run
#' estimation (see Examples below). The package deploys OpenMP that parallelizes computation of the Nadaraya-Watson
#' estimator. The default value of the number of threads is \code{parallel::detectCores()-1L}. To set
#' manually, please use \code{\link{set_num_threads}(x)}. If it's set to be 1L, the function doesn't use OpenMP.
#'
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
#' # identifiable set of parameters
#' ests_true <- c(1, -.5)
#'
#' # using matrix/vector
#' qmle0 <- semiBRM(x = X, y = Y, control = list(iterlim=50))
#'
#' # using formula and data
#' data <- data.frame(Y, X1, X2, X3)
#' qmle1 <- semiBRM(Y~X1+X2+X3, data = data, control = list(iterlim=50))
#'
#' @export
semiBRM <- function(x, y, ...) UseMethod("semiBRM")


#' @rdname semiBRM
#' @export
semiBRM.default <- function(x, y, r = 1/6.01, tau = 0.025, ...)
{
  if (is.factor(y)){
    y <- as.integer(levels(y))[y]
  }else if (is.numeric(y)){
    y <- as.integer(round(y))
  }else if (is.character(y)){
    y <- match(y, unique(y))-1L
  }

  if (sum(y %in% c(0L, 1L) == FALSE) > 0){
    stop("'y' is not a binary variable.")
  }

  x <- as.matrix(x)
  N <- length(y)

  # log likelihood function
  logLik <- function(parms)
  {
    v <- x[,1L] + as.vector(x[,-1L]%*%parms)
    h <- stats::sd(v)*1.06*N^(-r)

    prob <- GaussianNadarayaWatsonEstimator(v, y, h)

    if ( is.null(trimming_indicator) ){
      loglik <- log( c(prob[y==1L], 1-prob[y==0L]) )
      return( sum(loglik) )
    }

    prob_ <- prob[trimming_indicator == 1L]
    y_ <- y[trimming_indicator == 1L]

    loglik <- log( c(prob_[y_ == 1L], 1-prob_[y_ == 0L]) )
    return( sum(loglik) )
  }

  # starting values
  bols <- stats::coefficients(RcppArmadillo::fastLm(cbind(x, 1), y))
  start <- (bols[2L:(length(bols)-1L)]/bols[1L])

  # pilot estimation
  trimming_indicator <- semiBRM::TrimmingIndicator(x, tau)
  qmle_0 <- maxLik::maxLik(logLik, start = start, method = "BFGS", ...)

  # main estimation
  vhat <- x[,1L] + as.vector(x[,-1L]%*%qmle_0$estimate)
  trimming_indicator <- semiBRM::TrimmingIndicator(vhat, tau)
  qmle_1 <- maxLik::maxLik(logLik, start = qmle_0$estimate, method = "BFGS", ...)

  # return
  output <- list(estimate = qmle_1$estimate,
                 log.likelihood = qmle_1$maximum,
                 gradient = qmle_1$gradient,
                 hessian = qmle_1$hessian,
                 code = qmle_1$code,
                 message = qmle_1$message,
                 iter = unname(qmle_1$iterations),
                 control = qmle_1$control,
                 model = data.frame(y, x),
                 r = r,
                 trimming.level = tau)

  output$call <- match.call()
  class(output) <- "semiBRM"

  output
}


#' @rdname semiBRM
#' @export
semiBRM.formula <- function(formula, data, r = 1/6.01, tau = 0.025, ...)
{
  form <- stats::update(formula, . ~ . - 1)
  mf <- stats::model.frame(formula = form, data = data)
  x <- stats::model.matrix(attr(mf, "terms"), data = mf)
  y <- stats::model.response(mf)

  est <- semiBRM.default(x = x, y = y, r = r, tau = tau, ...)
  est$call <- match.call()
  est$formula <- formula
  est
}


#' Predict in-sample semiparametric conditional probability
#'
#' This calculates in-sample pointwise semiparametric conditional probability \code{Pr(y=1|x)} based
#' on the fitted semiparametric binary response model. It uses the leave-one-out version of
#' the Nadaraya-Watson estimator.
#'
#' @param object a fitted 'semiBRM' object.
#' @param ... further arguments (currently ignored).
#'
#' @return a numeric vector of in-sample pointwise semiparametric conditional probability.
#'
#' @rdname fitted
#' @export
fitted.semiBRM <- function(object, ...)
{
  vhat <- object$model[,2L] + as.vector(as.matrix(object$model[,-c(1L,2L)])%*%object$estimate)
  h <- stats::sd(vhat)*1.06*length(vhat)^(-1/5)

  GaussianNadarayaWatsonEstimator(vhat, object$model[,1L], h)
}

#' Predict semiparametric conditional probability
#'
#' This calculates pointwise semiparametric conditional probability \code{Pr(y=1|x)} based
#' on the fitted semiparametric binary response model. If new data are not given, it returns
#' in-sample conditional probabilities using \code{\link{fitted.semiBRM}}.
#'
#' @param object a fitted 'semiBRM' object.
#' @param newdata a data.frame or matrix at which conditional probabilities are computed.
#' @param boot.se a logical indicating whether to report standard errors and confidence intervals.
#' If \code{boot.se = TRUE}, it calculates standard errors pointwise from semiparametric bootstrapping.
#' @param ci.level a numeric representing bootstrap confidence intervals. This is useful only when
#' \code{boot.se = TRUE}.
#' @param nboot an integer indicating the number of bootstrap replications. This is useful only when
#' \code{boot.se = TRUE}.
#' @param ... further arguments (currently ignored).
#'
#' @return
#'  If \code{boot.se = FALSE}, then it will return a list with the following components:
#'
#' | prob | predictions.|
#' | ---- | ----------- |
#' | non.endpoint | taking \code{TRUE} if estimated probabilities are evaluated at points away from boundaries and \code{FALSE} otherwise. |
#'
#' If \code{boot.se = TRUE}, then it will return a list with the following components:
#'
#' | prob | predictions.|
#' | ---- | ----------- |
#' | boot.se | semiparametric bootstrap standard errors. |
#' | boot.ci | semiparametric bootstrap confidence intervals. |
#' | non.endpoint | taking \code{TRUE} if estimated probabilities are evaluated at points away from boundaries and \code{FALSE} otherwise. |
#'
#' @rdname predict
#' @export
predict.semiBRM <- function(object, newdata = NULL, boot.se = FALSE, ci.level = 0.95, nboot = 300L, ...)
{
  if (is.null(newdata)){

    semi_prob <- fitted.semiBRM(object)
    vhat <- object$model[,2L] + as.vector(as.matrix(object$model[,-c(1L,2L)])%*%object$estimate)
    trimming <- TrimmingIndicator(vhat, object$trimming.level)

    if (boot.se == TRUE){
      N <- length(vhat)
      h <- stats::sd(vhat)*1.06*N^(-1/5)
    }

  }else{

    if ( !is.null(object$formula) & !is.matrix(newdata) ){
      form <- stats::update(object$formula, NULL ~ . -1)
      args <- stats::model.matrix(form, newdata)
    }else{
      args <- newdata
    }

    stopifnot(is.matrix(args))

    vnew <- args[,1L] + as.vector(args[,-1L]%*%object$estimate)
    vhat <- object$model[,2L] + as.vector(as.matrix(object$model[,-c(1L,2L)])%*%object$estimate)

    N <- length(vhat)
    h <- stats::sd(vhat)*1.06*N^(-1/5)

    semi_prob <- GaussianNadarayaWatsonEstimator(vhat, object$model[,1L], h, vnew)

    bounds <- stats::quantile(vhat, c(object$trimming.level, 1-object$trimming.level))
    trimming <- ifelse(bounds[1L] < vnew & vnew < bounds[2L], 1L, 0L)
  }

  if (boot.se == TRUE){

    # turn off multithreading temporarily
    cur_num_threads <- get_num_threads()
    set_num_threads(1L)

    alpha <- (1-ci.level)/2
    Nargs <- length(semi_prob)

    boot_se <- rep(NaN, Nargs)
    boot_ci <- matrix(NaN, nrow = Nargs, ncol = 2L)
    colnames(boot_ci) <- c("lower", "upper")

    for (i in seq_len(Nargs)){

      y_boot <- matrix(ifelse(stats::runif(N*nboot) <= semi_prob[i], 1, 0), nrow = N, ncol = nboot)

      prob_boot <- rep(NaN, nboot)
      for (b in seq_len(nboot)){

        if (is.null(newdata)){
          prob_boot[b] <- GaussianNadarayaWatsonEstimator(vhat[-i], y_boot[-i,b], h, vhat[i])
        }else{
          prob_boot[b] <- GaussianNadarayaWatsonEstimator(vhat, y_boot[,b], h, vnew[i])
        }
      }

      boot_se[i] <-stats::sd(prob_boot)
      boot_ci[i,] <-stats::quantile(prob_boot, c(alpha, 1-alpha))
    }

    output <- list("prob" = semi_prob,
                   "boot.se" = boot_se,
                   "boot.ci" = boot_ci,
                   "non.endpoint" = as.logical(trimming))

    # restore thread setting
    set_num_threads(cur_num_threads)

  }else{
    output <- list("prob" = semi_prob,
                   "non.endpoint" = as.logical(trimming))

  }

  output
}


#' Calculate asymptotic variance-covariance matrix
#'
#' This calculates the asymptotic variance-covariance matrix of the estimates in the
#' fitted semiparametric binary response model as the inverse of the negative hessian matrix.
#'
#' @param object a fitted 'semiBRM' object.
#'
#' @param ... further arguments (currently ignored).
#'
#' @return a numeric matrix of the asymptotic variance-covariance.
#'
#' @rdname vcov
#' @export
vcov.semiBRM <- function(object, ...)
{
  -solve(object$hessian)
}


#' Extract parameter estimates
#'
#' This returns parameter estimates of the fitted semiparametric binary response model.
#'
#' @param object a fitted 'semiBRM' object.
#' @param ... further arguments (currently ignored).
#'
#' @return a numeric vector of parameter estimates.
#'
#' @rdname coef
#' @export
coef.semiBRM <- function(object, ...)
{
  object$estimate
}


#' @export
print.semiBRM <- function(x, ...)
{
  cat("Quasi maximum likelihood estimation using 'maxLik' with BFGS method", fill = TRUE)
  cat(sprintf("Return code %d: %s(%d iterations)", x$code, x$message, x$iter), fill = TRUE)
  cat(sprintf("Log-likelihood: %.4f (%d free parameter(s))", x$log.likelihood, length(x$estimate)), fill = TRUE)
  cat(sprintf("Estimates: %s", paste0(sprintf("%.6f", x$estimate), collapse = ", ")), fill = TRUE)
}


#' Summary of quasi maximum likelihood estimation
#'
#' This shows summary of quasi maximum likelihood estimation of the parameters in the fitted
#' semiparametric binary response model.
#'
#' @param object a fitted 'semiBRM' object.
#' @param ... further arguments (currently ignored).
#'
#' @rdname summary
#' @export
summary.semiBRM <- function(object, ...)
{
  vcov <- -solve(object$hessian)
  se <- sqrt(diag(vcov))
  tval <- object$estimate / se

  coefs <- cbind(Estimate = c(1, object$estimate),
                 Std.Error = c(NaN, se),
                 t.value = c(NaN, tval),
                 p.value = c(NaN, 2*stats::pnorm(-abs(tval))))

  rownames(coefs) <- names(object$model)[-1L]

  res <- list(call = object$call,
              coefficients = coefs,
              code = object$code,
              message = object$message,
              iter = object$iter,
              log.likelihood = object$log.likelihood,
              obs = nrow(object$model))

  class(res) <- "summary.semiBRM"
  res
}


#' @export
print.summary.semiBRM <- function(x, ...)
{
  cat("Semiparametric Binary Response Model", fill = TRUE)
  cat("\n")

  cat("Call:", fill = TRUE)
  print(x$call)

  cat("\nEstimates:", fill = TRUE)
  stats::printCoefmat(x$coefficients, P.value = TRUE, has.Pvalue = TRUE, na.print = "")
  cat(sprintf("\nLog-likelihood: %.4f (%d free parameter(s), Obs: %d)",
              x$log.likelihood, nrow(x$coefficients)-1L, x$obs), fill = TRUE)

  cat("\nNote:", fill = TRUE)
  cat("The coefficient of the first explanatory variable is normalized", fill = TRUE)
  cat("to 1, and the rest of the coefficients are rescaled conformably", fill = TRUE)
  cat("to it, i.e., the estimands are the rescaled coefficients.", fill = TRUE)
}

