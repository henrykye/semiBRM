# data generating process -----------------------------------------------------------------------------------------
library(semiBRM)

N <- 1500L
X1 <- rnorm(N)
X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
X3 <- rnorm(N)^2/sqrt(2)
X <- cbind(X1, X2, X3)
beta <- c(2, 2, -1, -1)
V <- as.vector(cbind(X, 1)%*%beta)
Y <- ifelse(V >= rnorm(N), 1L, 0L)


# parameter estimation --------------------------------------------------------------------------------------------

# true parameters
coefs_true <- c(1, -.5)

# semiparametric approach
qmle <- semiBRM(X = X, Y = Y, control = list(iterlim=50))
coefs_semi <- qmle$estimate

# parametric approach (probit)
probit <- glm(Y ~ X, family = binomial(link = "probit"))
coefs_probit <- probit$coefficients[-1L][-1L]/probit$coefficients[-1L][1L]


# conditional probability -----------------------------------------------------------------------------------------

# true conditional probability at the means
Xbar <- colMeans(X)
Prob_0 <- pnorm( as.vector(c(Xbar, 1)%*%beta) )

# semiparametric approach
Vhat <- X[,1L] + as.vector(X[,-1L]%*%semi_coefs)
Vhat_bar <- Xbar[1L] + as.vector(Xbar[-1L]%*%semi_coefs)
h <- sd(Vhat)*N^(-1/5)
Prob_semi <- GaussianKerNonpar(Y, Vhat, Vhat_bar, h)

# parametric approach
Prob_prob <- pnorm(as.vector(probit$coefficients%*%c(1, Xbar)))


# outputting ------------------------------------------------------------------------------------------------------
{
cat(sprintf("        %7s\t%7s\t%7s", "true", "Semi", "Probit"), fill = TRUE)
cat(sprintf("parm 1: %7.4f\t%7.4f\t%7.4f", coefs_true[1L], coefs_semi[1L], coefs_probit[1L]), fill = TRUE)
cat(sprintf("parm 2: %7.4f\t%7.4f\t%7.4f", coefs_true[2L], coefs_semi[2L], coefs_probit[2L]), fill = TRUE)
cat(sprintf("prob. : %7.4f\t%7.4f\t%7.4f", Prob_0, Prob_semi, Prob_prob), fill = TRUE)
}

