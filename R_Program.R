############################################################
## FULL SIMULATION CODE FOR SEMIPARAMETRIC SIEVE ESTIMATION
############################################################

############################
# Required Packages
############################
library(splines)
library(MASS)
library(Matrix)
library(numDeriv)
library(stats)

############################
# True Model Specification
############################

beta.true <- 1.5

m.true <- function(z){
  sin(2*pi*z)
}

############################
# Data Generation
############################

generate.data <- function(n, regime = "mild"){
  
  X <- rnorm(n)
  Z <- runif(n,0,1)
  
  mu <- beta.true*X + m.true(Z)
  
  eps <- rnorm(n,0,0.5)
  
  Y.full <- mu + eps
  
  ## Truncation variable
  if(regime=="mild"){
    T <- rnorm(n, mean = -1, sd = 1)
  }
  
  if(regime=="moderate"){
    T <- rnorm(n, mean = 0, sd = 1)
  }
  
  if(regime=="severe"){
    T <- rnorm(n, mean = 1, sd = 1)
  }
  
  delta <- as.numeric(Y.full >= T)
  
  Y.obs <- ifelse(delta==1, Y.full, NA)
  
  list(
    X=X,
    Z=Z,
    Y=Y.obs,
    delta=delta,
    Y.full=Y.full
  )
}

############################
# B-spline Basis
############################

build.basis <- function(z, K){
  
  bs(
    z,
    df = K,
    degree = 3,
    intercept = TRUE
  )
}

############################
# Penalized Estimation
############################

fit.semiparametric <- function(X,Z,Y,delta,K,lambda){
  
  obs <- which(delta==1)
  
  Xo <- X[obs]
  Zo <- Z[obs]
  Yo <- Y[obs]
  
  B <- build.basis(Zo,K)
  
  D <- diff(diag(K), differences=2)
  
  P <- t(D)%*%D
  
  XtX <- sum(Xo^2)
  XtB <- t(Xo)%*%B
  BtB <- t(B)%*%B + lambda*P
  
  A11 <- matrix(XtX,1,1)
  A12 <- XtB
  A21 <- t(A12)
  A22 <- BtB
  
  A <- rbind(
    cbind(A11,A12),
    cbind(A21,A22)
  )
  
  rhs <- c(
    sum(Xo*Yo),
    t(B)%*%Yo
  )
  
  conv <- TRUE
  
  est <- tryCatch(
    solve(A,rhs),
    error = function(e){
      conv <<- FALSE
      rep(NA,K+1)
    }
  )
  
  beta.hat <- est[1]
  gamma.hat <- est[-1]
  
  list(
    beta = beta.hat,
    gamma = gamma.hat,
    conv = conv,
    basis = B
  )
}

############################
# L2 Error Calculation
############################

L2.error <- function(gamma,Z,K){
  
  grid <- seq(0,1,length=300)
  
  B.grid <- build.basis(grid,K)
  
  m.hat <- B.grid %*% gamma
  
  truth <- m.true(grid)
  
  sqrt(
    mean(
      (m.hat-truth)^2
    )
  )
}

############################
# Confidence Interval
############################

compute.se <- function(X,Z,delta,K,lambda){
  
  obs <- which(delta==1)
  
  Xo <- X[obs]
  Zo <- Z[obs]
  
  B <- build.basis(Zo,K)
  
  D <- diff(diag(K), differences=2)
  
  P <- t(D)%*%D
  
  XtX <- sum(Xo^2)
  XtB <- t(Xo)%*%B
  BtB <- t(B)%*%B + lambda*P
  
  A11 <- matrix(XtX,1,1)
  A12 <- XtB
  A21 <- t(A12)
  A22 <- BtB
  
  A <- rbind(
    cbind(A11,A12),
    cbind(A21,A22)
  )
  
  vcov.mat <- tryCatch(
    solve(A),
    error=function(e) NULL
  )
  
  if(is.null(vcov.mat)){
    return(NA)
  }
  
  sqrt(vcov.mat[1,1])
}

############################
# Single Simulation Run
############################

single.run <- function(
    n,
    regime,
    K,
    lambda
){
  
  dat <- generate.data(n, regime)
  
  fit <- fit.semiparametric(
    dat$X,
    dat$Z,
    dat$Y,
    dat$delta,
    K,
    lambda
  )
  
  if(!fit$conv){
    return(NULL)
  }
  
  se.beta <- compute.se(
    dat$X,
    dat$Z,
    dat$delta,
    K,
    lambda
  )
  
  cover <- as.numeric(
    (fit$beta - 1.96*se.beta <= beta.true) &
      (fit$beta + 1.96*se.beta >= beta.true)
  )
  
  l2 <- L2.error(
    fit$gamma,
    dat$Z,
    K
  )
  
  list(
    beta = fit$beta,
    l2 = l2,
    cover = cover
  )
}

############################
# Monte Carlo Simulation
############################

run.simulation <- function(
    n,
    regime,
    M = 500,
    K.type = "optimal",
    lambda.type = "optimal"
){
  
  ## Sieve dimension
  if(K.type=="optimal"){
    K <- floor(n^(1/5))
  }
  
  if(K.type=="large"){
    K <- floor(n^0.6)
  }
  
  ## Penalty
  if(lambda.type=="optimal"){
    lambda <- n^(-1)
  }
  
  if(lambda.type=="small"){
    lambda <- n^(-3)
  }
  
  beta.vec <- c()
  l2.vec <- c()
  cover.vec <- c()
  
  conv.fail <- 0
  
  for(m in 1:M){
    
    out <- single.run(
      n=n,
      regime=regime,
      K=K,
      lambda=lambda
    )
    
    if(is.null(out)){
      conv.fail <- conv.fail + 1
      next
    }
    
    beta.vec <- c(beta.vec,out$beta)
    l2.vec <- c(l2.vec,out$l2)
    cover.vec <- c(cover.vec,out$cover)
  }
  
  bias <- mean(beta.vec-beta.true)
  sdv <- sd(beta.vec)
  rmse <- sqrt(mean((beta.vec-beta.true)^2))
  coverage <- mean(cover.vec)*100
  l2.avg <- mean(l2.vec)
  
  list(
    Bias = bias,
    SD = sdv,
    RMSE = rmse,
    Coverage = coverage,
    L2 = l2.avg,
    Convergence = 100*(1-conv.fail/M)
  )
}

############################
# Main Simulation Design
############################

n.vec <- c(200,500,1000)

regimes <- c(
  "mild",
  "moderate",
  "severe"
)

results <- list()

############################
# Main Results
############################

for(n in n.vec){
  
  for(reg in regimes){
    
    cat("\nRunning:", n, reg,"\n")
    
    results[[paste(n,reg,sep="_")]] <-
      run.simulation(
        n=n,
        regime=reg,
        M=500,
        K.type="optimal",
        lambda.type="optimal"
      )
  }
}

############################
# Effect of Sieve Dimension
############################

sieve.optimal <- run.simulation(
  n=500,
  regime="mild",
  M=500,
  K.type="optimal",
  lambda.type="optimal"
)

sieve.large <- run.simulation(
  n=500,
  regime="mild",
  M=500,
  K.type="large",
  lambda.type="optimal"
)

############################
# Effect of Undersmoothing
############################

lambda.optimal <- run.simulation(
  n=500,
  regime="mild",
  M=500,
  K.type="optimal",
  lambda.type="optimal"
)

lambda.small <- run.simulation(
  n=500,
  regime="mild",
  M=500,
  K.type="optimal",
  lambda.type="small"
)

############################
# Print Results
############################

cat("\n========================\n")
cat("MAIN RESULTS\n")
cat("========================\n")
print(results)

cat("\n========================\n")
cat("SIEVE EFFECT\n")
cat("========================\n")
print(sieve.optimal)
print(sieve.large)

cat("\n========================\n")
cat("UNDERSMOOTHING EFFECT\n")
cat("========================\n")
print(lambda.optimal)
print(lambda.small)