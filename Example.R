# setup ---------------------------------------------------------------------------------------
library(semiBRM)
set.seed(20190815)

# data generating process ----------------------------------------------------------------------

N <- 1500L
X1 <- rnorm(N)
X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
X3 <- rnorm(N)^2/sqrt(2)
X <- cbind(X1, X2, X3)
beta <- c(2, 2, -1, -1)
V <- as.vector(cbind(X, 1)%*%beta)
Y <- ifelse(V >= rnorm(N), 1L, 0L)


# parameter estimation -------------------------------------------------------------------------

## estimands: the rescaled coefficients by the first coefficient excluding intercept
coefs_true <- c(1, -.5)
data <- data.frame(Y, X1, X2, X3)

## Klein and Spady (1993): semiparametric approach
semi <- semiBRM(Y~X1+X2+X3, data = data, control = list(iterlim=50))
coefs_semi <- coef(semi)

## Probit: parametric approach
probit <- glm(Y ~ X1+X2+X3, family = binomial(link = "probit"), data = data)
coefs_probit <- probit$coefficients[-1L][-1L]/probit$coefficients[-1L][1L]

## formatted print
{
    cat(sprintf("        %7s %7s %7s\n", "True", "Probit", "Semi"))
    cat(sprintf("parm 1: %7.4f %7.4f %7.4f\n", coefs_true[1L], coefs_probit[1L], coefs_semi[1L]))
    cat(sprintf("parm 2: %7.4f %7.4f %7.4f\n", coefs_true[2L], coefs_probit[2L], coefs_semi[2L]))
}



# in-sample prediction ------------------------------------------------------------------------

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


# out-of-sample prediction --------------------------------------------------------------------

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

# average marginal effects --------------------------------------------------------------------

## marginal effects as difference betwen conditional probabilities with and without perturbation
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


# quantitle marginal effects ------------------------------------------------------------------

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

