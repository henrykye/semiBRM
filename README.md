
<!-- README.md is generated from README.Rmd. Please edit that file -->

# semiBRM: Semiparametric Binary Response Models in R

This **R** package **semiBRM** offers an implementation of single-index
semiparametric binary response models, theorized in a seminal work in
the semiparametric econometrics literature [Klein and Spady
(1993)](https://doi.org/10.2307/2951556).

This is built upon [Rcpp](http://www.rcpp.org) along with OpenMP, which
parallelizes computation of nonparametric conditional expectations over
given data points. This improves computation efficiency enormously in
parameter estimation, which is one of the main features of this package.

The package is still in development, with the goal of making it as easy
to employ [Klein and Spady (1993)](https://doi.org/10.2307/2951556) as
to run Probit or Logistic regressions in **R**.

## Identifiable Set of Parameters

In a single-index semiparametric approach, not all parameters are
identifiable. First of all, intercept cannot be estimated. Secondly, the
coefficients can be estimated as ratios to a “basis” coefficient. This
package sets the coefficient of the first explanatory variable as the
basis so that its coefficient is normalized to 1, and those of the rest
of the variables are estimated conformably to it.

In the context of binary response models, or nonlinear models more
broadly, this is not problematic at all, as the coefficients themselves
have little interpretation. As the matter of fact, the identifiable set
of parameters in this approach can correctly estimate the conditional
probability of Pr(Y=1|X), which would be the quantity of interest in
many cases.

## Installation

The development version can be installed from GitHub:

``` r
library(devtools)
install_github(repo="henrykye/semiBRM")
```

### This version was built with:

  - R version 4.0.5 (2021-03-31)
  - maxLik 1.4-8
  - Rcpp 1.0.6
  - RcppArmadillo 0.10.4.0.0

## Example

### Setup

``` r
library(semiBRM)
set.seed(20190815) # for reproduction of results
```

### Data Generating Process

``` r
## data generating process
N <- 1500L
X1 <- rnorm(N)
X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
X3 <- rnorm(N)^2/sqrt(2)
X <- cbind(X1, X2, X3)
beta <- c(2, 2, -1, -1) # this is the original set of coefficients
V <- as.vector(cbind(X, 1)%*%beta)
Y <- ifelse(V >= rnorm(N), 1L, 0L)
```

### Parameter Estimation

``` r
## estimands: the rescaled coefficients by the first coefficient excluding intercept
coefs_true <- c(1, -.5)
data <- data.frame(Y, X1, X2, X3)

## Klein and Spady (1993): semiparametric approach
semi <- semiBRM(Y ~ X1 + X2 + X3, data = data, control = list(iterlim=50))
coefs_semi <- coef(semi)

## Probit: parametric approach
probit <- glm(Y ~ X1 + X2 + X3, family = binomial(link = "probit"), data = data)
coefs_probit <- probit$coefficients[-1L][-1L]/probit$coefficients[-1L][1L]

## formatted print
{
    cat(sprintf("        %7s %7s %7s\n", "True", "Probit", "Semi"))
    cat(sprintf("parm 1: %7.4f %7.4f %7.4f\n", coefs_true[1L], coefs_probit[1L], coefs_semi[1L]))
    cat(sprintf("parm 2: %7.4f %7.4f %7.4f\n", coefs_true[2L], coefs_probit[2L], coefs_semi[2L]))
}
#>            True  Probit    Semi
#> parm 1:  1.0000  0.9296  0.9430
#> parm 2: -0.5000 -0.4069 -0.4103
```

### In-Sample Prediction

``` r
## in-sample conditional probability
in_prob_true <- pnorm(V)
in_prob_semi <- predict(semi)
in_prob_probit <- fitted(probit)

## formatted print
target <- sample.int(N, size = 10L)
{
    cat(sprintf("%7s %8s %8s %8s %12s\n", "Obs.", "True", "Probit", "Semi", "non.endpoint") )
    for (i in target){
        cat(sprintf("[%04d]: %8.6f %8.6f %8.6f %12s\n",
                    i, in_prob_true[i], in_prob_probit[i], in_prob_semi$prob[i], in_prob_semi$non.endpoint[i]))
    }
}
#>    Obs.     True   Probit     Semi non.endpoint
#> [0034]: 0.999400 0.999262 0.991678         TRUE
#> [0514]: 0.898871 0.859331 0.774815         TRUE
#> [0480]: 0.995355 0.994697 0.970159         TRUE
#> [1285]: 0.004911 0.006525 0.034459         TRUE
#> [0065]: 0.990712 0.993818 0.966249         TRUE
#> [0966]: 0.896823 0.862113 0.780052         TRUE
#> [0558]: 0.997327 0.997560 0.981306         TRUE
#> [1001]: 0.000001 0.000000 0.000653         TRUE
#> [0843]: 0.998430 0.998133 0.984302         TRUE
#> [0381]: 0.781482 0.773341 0.697235         TRUE
```

### Out-of-Sample Prediction

``` r
## conditional probability at the means
Xbar <- colMeans(X)
newdata <- as.data.frame(as.list(Xbar))

## predictions
out_prob_true <- pnorm(as.vector(c(Xbar, 1)%*%beta))
out_prob_semi <- predict(semi, newdata, boot.se = TRUE)
out_prob_probit <- pnorm(as.vector(coef(probit)%*%c(1, Xbar)))

## standard errors of Probit
grad <- dnorm(as.vector(coef(probit)%*%c(1, Xbar))) * c(1, Xbar)
out_stde_probit <- sqrt(as.vector(crossprod(grad, vcov(probit))%*%grad))

## formatted print
{
    cat(sprintf("            %7s %7s %7s\n", "True", "Probit", "Semi"))
    cat(sprintf("Prob. Est.: %7.4f %7.4f %7.4f\n", out_prob_true, out_prob_probit, out_prob_semi$prob))
    cat(sprintf(" Std. Err.: %7s %7.4f %7.4f\n", "", out_stde_probit, out_prob_semi$boot.se))
}
#>                True  Probit    Semi
#> Prob. Est.:  0.5833  0.5792  0.5543
#>  Std. Err.:          0.0230  0.0221
```

### Average Marginal Effects

``` r
## marginal effects as difference between conditional probabilities with and without perturbation
delta <- sd(X1) # size of perturbation

me_true <- pnorm(as.vector(cbind(X1+delta, X2, X3, 1)%*%beta)) - pnorm(as.vector(cbind(X, 1)%*%beta))

## average marginal effects
av_me_true <- mean(me_true)
av_me_semi <- MarginalEffects(semi, variable = "X1", delta = delta)

## formatted print
{
    cat(sprintf("           %7s %7s \n", "True", "Semi"))
    cat(sprintf("Marg.Eff.: %7.4f %7.4f\n", av_me_true, av_me_semi[1L]))
    cat(sprintf("Std. Err.: %7s %7.4f\n", "", av_me_semi[2L]))
}
#>               True    Semi 
#> Marg.Eff.:  0.2021  0.2029
#> Std. Err.:          0.0050
```

### Quantile Marginal Effects

``` r
## percentile cutoffs
p.cutoffs <- c(1/4, 2/4, 3/4)

## group indicators
q_vals <- quantile(X1, p.cutoffs)
q1 <- ifelse(X1 <= q_vals[1L], 1L, 0L)
q2 <- ifelse(q_vals[1L] < X1 & X1 <= q_vals[2L], 1L, 0L)
q3 <- ifelse(q_vals[2L] < X1 & X1 <= q_vals[3L], 1L, 0L)
q4 <- ifelse(q_vals[3L] < X1, 1L, 0L)

## quantile marginal effects
q_me_true <- c("G1" = mean(me_true[q1==1]), "G2" = mean(me_true[q2==1]),
               "G3" = mean(me_true[q3==1]), "G4" = mean(me_true[q4==1]))

q_me_semi <- MarginalEffects(semi, variable = "X1", delta, p.cutoffs)

## formatted print
{
    cat("         ", sprintf(" %5s  ", paste0("G", 1:4)), fill = TRUE)
    cat("True    :", sprintf(" %.4f ", q_me_true), fill = TRUE)
    cat("Semi    :", sprintf(" %.4f ", q_me_semi[,1L]), fill = TRUE)
    cat("Std.Err.:", sprintf("(%.4f)", q_me_semi[,2L]), fill = TRUE)
}
#>               G1       G2       G3       G4  
#> True    :  0.2232   0.3175   0.2196   0.0481 
#> Semi    :  0.2186   0.3130   0.2248   0.0552 
#> Std.Err.: (0.0100) (0.0088) (0.0096) (0.0061)
```

## References

Eddelbuettel, D., & François, R. (2011). Rcpp: Seamless R and C++
integration. *Journal of Statistical Software*, *40(8)*, 1-18.
*<https://dirk.eddelbuettel.com/code/rcpp/Rcpp-introduction.pdf>*

Klein, R. W., & Spady, R. H. (1993). An Efficient Semiparametric
Estimator for Binary Response Models. *Econometrica*, *61(2)*, 387-421.
*<https://doi.org/10.2307/2951556>*.

Klein, R., & Vella, F. (2009). A Semiparametric Model for Binary
Response and Continuous Outcomes Under Index Heteroscedasticity.
*Journal of Applied Econometrics*, *24(5)*, 735-762.
*<https://doi.org/10.1002/jae.1064>*
