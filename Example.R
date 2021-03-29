
library(semiBRM)
# load_all()

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

# estimands: the rescaled coefficients by the first coefficient excluding intercept
coefs_true <- c(1, -.5)
data <- data.frame(Y, X1, X2, X3)

# Klein and Spady (1993): semiparametric approach
semi <- semiBRM(Y~X1+X2+X3, data = data, control = list(iterlim=50))
coefs_semi <- coef(semi)

# Probit: parametric approach
probit <- glm(Y ~ X1+X2+X3, family = binomial(link = "probit"), data = data)
coefs_probit <- probit$coefficients[-1L][-1L]/probit$coefficients[-1L][1L]

cbind(true=coefs_true, probit=coefs_probit, semi=coefs_semi)


# conditional probability ----------------------------------------------------------------------

# in-sample prediction of conditional probability
in_prob_true <- pnorm(V)
in_prob_semi <- predict(semi)
in_prob_probit <- fitted(probit)

round(cbind(true=in_prob_true, probit=in_prob_probit,
            semi=in_prob_semi$prob, non.endpoint=in_prob_semi$non.endpoint)[sample.int(N, 20),], 6)

# out-of-sample prediction of conditional probability at the means
Xbar <- colMeans(X)
newdata <- as.data.frame(as.list(Xbar))

out_prob_true <- pnorm(as.vector(c(Xbar, 1)%*%beta))
out_prob_semi <- predict(semi, newdata)
out_prob_probit <- pnorm(as.vector(probit$coefficients%*%c(1, Xbar)))

c(true=out_prob_true, probit=out_prob_probit, semi=out_prob_semi$prob)



# average marginal effect ---------------------------------------------------------------------


as.vector(cbind(1, X1, X2, X3)%*%probit$coefficients)

av_me_true <- mean(
    pnorm(as.vector(cbind(X1+sd(X1), X2, X3, 1)%*%beta)) - pnorm(as.vector(cbind(X, 1)%*%beta))
)

av_me_semi <-  MarginalEffects(semi, variable = "X1", delta = sd(X1))[1,1]

av_me_probit <- mean(
    pnorm(as.vector(cbind(1, X1+sd(X1), X2, X3)%*%coef(probit))) - pnorm(as.vector(cbind(1, X)%*%coef(probit)))
)

Vhat <- X[,1] + as.vector(X[,-1L]%*%coef(semi))
h <- sd(Vhat) * N^(-1/5)

Vnew <- X[,1] + sd(X[,1]) + as.vector(X[,-1L]%*%coef(semi))


mean(GaussianNadarayaWatsonEstimator(Vhat, Y, h, Vnew)-
         GaussianNadarayaWatsonEstimator(Vhat, Y, h))




# outputting ----------------------------------------------------------------------------------
{
cat(sprintf("        %7s\t%7s\t%7s", "true", "Semi", "Probit"), fill = TRUE)
cat(sprintf("parm 1: %7.4f\t%7.4f\t%7.4f", coefs_true[1L], coefs_semi[1L], coefs_probit[1L]), fill = TRUE)
cat(sprintf("parm 2: %7.4f\t%7.4f\t%7.4f", coefs_true[2L], coefs_semi[2L], coefs_probit[2L]), fill = TRUE)
cat(sprintf("prob. : %7.4f\t%7.4f\t%7.4f", Prob_0, Prob_semi, Prob_prob), fill = TRUE)
}

