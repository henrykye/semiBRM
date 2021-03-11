
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semiBRM: Semiparametric Binary Response Model

This R package `semiBRM` offers an implementation of single-index
semiparametric binary response models introduced in [Klein and Spady
(1993)](https://doi.org/10.2307/2951556).

This is built upon `Rcpp` along with `OpenMP`, which parallelizes
computation of nonparametric conditional expectations over given data
points. This helps improve computation efficiency enormously in
parameter estimation, which is one of the main features of this package.

This is still in development. The ultimate goal of it is to help
researchers employ [Klein and Spady
(1993)](https://doi.org/10.2307/2951556) as easy as running Probit or
Logistic regressions in R.

## Identifiable Parameters in Semiparametric Models

In a single-index semiparametric approach, not all parameters are
identifiable. First of all, intercept cannot be estimated. Secondly, the
coefficients can be estimated as ratios to a “basis” coefficient. This
package sets the coefficient of the first explanatory variable as the
basis so that its coefficient is normalized to 1, and those of the rest
of the variables are estimated conformably to it.

In the context of binary response models, or nonlinear models more
broadly, this is not problematic at all, as the coefficients themselves
have little meaning. As the matter of fact, the identifiable set of
parameters in this approach can correctly estimate the conditional
probability of Y being 1 given X, which would be the ultimate quantity
of interest in most of the cases.

## Installation

The development version can be installed from GitHub as follows:

``` r
library(devtools)
install_github(repo="hk599/semiBRM")
```

### This version was built with:

  - R 4.0.4
  - rtools40
  - maxLik 1.4-6
  - matrixStats 0.58.0
  - Rcpp 1.0.6
  - RcppArmadillo 0.10.2.2.0

## Example

### Implementation

``` r
library(semiBRM)

# data generating process
N <- 1500L
X1 <- rnorm(N)
X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
X3 <- rnorm(N)^2/sqrt(2)
X <- cbind(X1, X2, X3)
beta <- c(2, 2, -1, -1) # this is an original set of coefficients
V <- as.vector(cbind(X, 1)%*%beta)
Y <- ifelse(V >= rnorm(N), 1L, 0L)


# parameter estimation
coefs_true <- c(1, -.5) # this is an identifiable set of coefficients

## semiparametric approach
qmle <- semiBRM(X = X, Y = Y, control = list(iterlim=50))
coefs_semi <- qmle$estimate

## parametric approach (probit)
probit <- glm(Y ~ X, family = binomial(link = "probit"))
coefs_probit <- probit$coefficients[-1L][-1L]/probit$coefficients[-1L][1L]


# conditional probability at the means
Xbar <- colMeans(X)
Prob_0 <- pnorm( as.vector(c(Xbar, 1)%*%beta) ) # true conditional probability at the means

## semiparametric approach
Vhat <- X[,1L] + as.vector(X[,-1L]%*%coefs_semi)
Vhat_bar <- Xbar[1L] + as.vector(Xbar[-1L]%*%coefs_semi)
h <- sd(Vhat)*N^(-1/5)
Prob_semi <- GaussianKerNonpar(Y, Vhat, Vhat_bar, h)

## parametric approach
Prob_prob <- pnorm(as.vector(probit$coefficients%*%c(1, Xbar)))
```

### Results

    #>            true     Semi  Probit
    #> parm 1:  1.0000   0.9536  0.9712
    #> parm 2: -0.5000  -0.4938 -0.5059
    #> prob. :  0.6217   0.6155  0.6478

## References

Klein, R. W., & Spady, R. H. (1993). An Efficient Semiparametric
Estimator for Binary Response Models. *Econometrica*, 61(2), 387-421.
*<https://doi.org/10.2307/2951556>*.
