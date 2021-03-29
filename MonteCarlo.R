library(semiBRM)

hhmmss <- function(seconds){
    hh = seconds%/%3600
    mm = (seconds - hh*3600)%/%60
    ss = seconds - hh*3600 - mm*60
    return(sprintf("%02.0f:%02.0f:%02.0f", hh, mm, ss))
}


# marginal effects ----------------------------------------------------------------------------

S <- 300L

N <- 1000L
beta <- c(2, 2, -1, -1)

coef <- matrix(NaN, nrow = S, 2L)
coef.stder <- matrix(NaN, nrow = S, 2L)

avme.true <- matrix(NaN, nrow = S, 1L)
avme.semi <- matrix(NaN, nrow = S, 1L)
avme.stder <- matrix(NaN, nrow = S, 1L)
avme.error <- matrix(NaN, nrow = S, 1L)

qme.true <- matrix(NaN, nrow = S, 3L)
qme.semi <- matrix(NaN, nrow = S, 3L)
qme.stder <- matrix(NaN, nrow = S, 3L)
qme.error <- matrix(NaN, nrow = S, 3L)

start <- proc.time()

for (s in seq_len(S)){

    X1 <- rnorm(N)
    X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
    X3 <- rnorm(N)^2/sqrt(2)
    X <- cbind(X1, X2, X3)
    V <- as.vector(cbind(X, 1)%*%beta)
    Y <- ifelse(V >= rnorm(N), 1L, 0L)

    qmle <- semiBRM(x = X, y = Y, control = list(iterlim=50))

    # coefficients
    coef[s,] <- coef(qmle)
    coef.stder[s,] <- sqrt(diag(vcov(qmle)))


    # average marginal effects
    me <- pnorm(as.vector(cbind(X1+sd(X1), X2, X3, 1)%*%beta)) -
        pnorm(as.vector(cbind(X, 1)%*%beta))

    av_me_true <- mean(me)

    av_me_semi <-  MarginalEffects(qmle, variable = "X1", delta = sd(X1))

    avme.true[s,] <- av_me_true
    avme.semi[s,] <- av_me_semi[1L]
    avme.stder[s,] <- av_me_semi[2L]
    avme.error[s,] <- av_me_semi[1L] - av_me_true

    # quantile marginal effects
    x1.cutoff <- quantile(X1, c(1/3, 2/3))

    g1_true <- mean(me[ X1 <= x1.cutoff[1L]])
    g2_true <- mean(me[ X1 >  x1.cutoff[1L] & X1 <= x1.cutoff[2L] ])
    g3_true <- mean(me[ X1 >  x1.cutoff[2L] ])

    q_me_semi <-  MarginalEffects(qmle, variable = "X1", delta = sd(X1), p.cutoffs = c(1/3, 2/3))

    qme.true[s,] <-  c(g1_true, g2_true, g3_true)
    qme.semi[s,] <- q_me_semi[,1L]
    qme.stder[s,] <- q_me_semi[,2L]
    qme.error[s,] <- q_me_semi[,1L] - c(g1_true, g2_true, g3_true)

    cat(sprintf("progress = %03d/%03d done (elapsed = %s)\r", s, S, hhmmss( (proc.time()-start)[3L] )))

    if (s==S) cat("\n");
}


# without trimming ----------------------------------------------------------------------------

# matrixStats::colMeans2(avme.true)
# [1] 0.2005231
# matrixStats::colMeans2(avme.semi)
# [1] 0.1940123
# matrixStats::colMeans2(avme.stder)
# [1] 0.006110306
# matrixStats::colSds(avme.error)
# [1] 0.00696816
# matrixStats::colMeans2(avme.error)
# [1] -0.006510768

# matrixStats::colMeans2(qme.true)
# [1] 0.24908598 0.27697798 0.07535945
# matrixStats::colMeans2(qme.semi)
# [1] 0.23737675 0.26203671 0.08249326
# matrixStats::colMeans2(qme.stder)
# [1] 0.010624265 0.010250780 0.007762536
# matrixStats::colSds(qme.error)
# [1] 0.012519107 0.011962739 0.004470873
# matrixStats::colMeans2(qme.error)
# [1] -0.011709231 -0.014941269  0.007133807



# point estimation: semiparametric bootstrapping ----------------------------------------------

S <- 300L

N <- 1200L
beta <- c(2, 2, -1, -1)
Xbar <- t(c(0, 1, 1/sqrt(2)))

probs_semi <- rep(NaN, S)
probs_boot <- rep(NaN, S)
stder_boot <- rep(NaN, S)

start <- proc.time()
for (s in seq_len(S)){

    X1 <- rnorm(N)
    X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
    X3 <- rnorm(N)^2/sqrt(2)
    X <- cbind(X1, X2, X3)
    V <- as.vector(cbind(X, 1)%*%beta)
    Y <- ifelse(V >= rnorm(N), 1L, 0L)

    qmle <- semiBRM(x = X, y = Y, control = list(iterlim=50))

    prob.true <- pnorm(as.vector(cbind(Xbar, 1)%*%beta))
    prob.semi <- predict(qmle, newdata = Xbar)

    Vhat <- X[,1L] + as.vector(X[,-1L]%*%coef(qmle))
    Vnew <- Xbar[,1L] + as.vector(Xbar[,-1L]%*%coef(qmle))
    h <- sd(Vhat)*N^(-1/5)

    ## semiparametric bootstrapping
    B <- 500L
    prob_boot <- rep(NaN, B)
    for (b in seq_len(B)){
        Y.boot <- ifelse(runif(N) <= prob.semi$prob, 1L, 0L)
        prob_boot[b] <- GaussianNadarayaWatsonEstimator(Vhat, Y.boot, h, Vnew)
    }

    probs_semi[s] <- prob.semi$prob
    probs_boot[s] <- mean(prob_boot)
    stder_boot[s] <- sd(prob_boot)

    cat(sprintf("progress = %03d/%03d done (elapsed = %s)\r", s, S, hhmmss( (proc.time()-start)[3L] )))
    if (s == S) cat("\n");
}

prob.true
mean(probs_semi)
mean(probs_boot)

sd(probs_semi)
mean(stder_boot)



# point estimation: parameter simulation ------------------------------------------------------

S <- 300L

N <- 1200L
beta <- c(2, 2, -1, -1)

Xbar <- t(c(0, 1, 1/sqrt(2)))

probs_semi <- rep(NaN, S)
error_semi <- rep(NaN, S)
stder_semi <- rep(NaN, S)


start <- proc.time()
for (s in seq_len(S)){

    X1 <- rnorm(N)
    X2 <- (X1 + 2*rnorm(N))/sqrt(5) + 1
    X3 <- rnorm(N)^2/sqrt(2)
    X <- cbind(X1, X2, X3)
    V <- as.vector(cbind(X, 1)%*%beta)
    Y <- ifelse(V >= rnorm(N), 1L, 0L)

    qmle <- semiBRM(x = X, y = Y, control = list(iterlim=50))

    prob.true <- pnorm(as.vector(cbind(Xbar, 1)%*%beta))
    prob.semi <- predict(qmle, newdata = Xbar)

    Vhat <- X[,1L] + as.vector(X[,-1L]%*%coef(qmle))
    h <- sd(Vhat)*N^(-1/5)

    parms <- cbind(1, mvrnorm(n = 100, mu=coef(qmle), Sigma=vcov(qmle)))
    Vnew <- as.vector(tcrossprod(parms, Xbar))

    prob.simul <- GaussianNadarayaWatsonEstimator(Vhat, Y, h, Vnew)

    probs_semi[s] <- prob.semi$prob
    stder_semi[s] <- sd(prob.simul)
    error_semi[s] <- prob.semi$prob - prob.true


    cat(sprintf("progress = %03d/%03d done (elapsed = %s)\r", s, S, hhmmss( (proc.time()-start)[3L] )))
    if (s==S) cat("\n");
}

mean(probs_semi)
mean(stder_semi)
sd(probs_semi)

