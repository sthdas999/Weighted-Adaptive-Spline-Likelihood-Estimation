# ============================================================
# Comparison:
# Proposed Truncated vs Conventional Untruncated Sieve Estimator
# ============================================================

rm(list = ls())

library(splines)

set.seed(12345)

# ------------------------------------------------------------
# 1. Simulation settings
# ------------------------------------------------------------

n.values <- c(200, 500, 1000)

R <- 1000

r <- 2

sigma <- 0.5

beta.true <- 1

# Truncation regimes
truncation.regimes <- list(
  Mild     = c(-2, 2),
  Moderate = c(-1.5, 1.5),
  Severe   = c(-1, 1)
)

# Number of test observations
N.test <- 5000

# ------------------------------------------------------------
# 2. True nonlinear function
# ------------------------------------------------------------

m.true <- function(z) {
  sin(2 * pi * z)
}

# ------------------------------------------------------------
# 3. Generate complete data
# ------------------------------------------------------------

generate.complete.data <- function(n) {
  
  Z <- runif(n, 0, 1)
  
  X <- rnorm(n, 0, 1)
  
  eps <- rnorm(n, 0, sigma)
  
  Y <- beta.true * X + m.true(Z) + eps
  
  data.frame(
    Y = Y,
    X = X,
    Z = Z
  )
}

# ------------------------------------------------------------
# 4. Apply response truncation
# ------------------------------------------------------------

apply.truncation <- function(dat, interval) {
  
  dat[
    dat$Y >= interval[1] &
      dat$Y <= interval[2],
  ]
  
}

# ------------------------------------------------------------
# 5. Spline dimension
# ------------------------------------------------------------

get.K <- function(n) {
  
  floor(n^(1 / (2 * r + 1)))
  
}

# ------------------------------------------------------------
# 6. Construct spline basis
# ------------------------------------------------------------

make.spline.basis <- function(Z, K) {
  
  # Internal knots
  if (K > 1) {
    
    knots <- seq(
      min(Z),
      max(Z),
      length.out = K + 1
    )[-c(1, K + 1)]
    
  } else {
    
    knots <- NULL
    
  }
  
  bs(
    Z,
    knots = knots,
    degree = 3,
    intercept = TRUE
  )
  
}

# ------------------------------------------------------------
# 7. Fit proposed truncated likelihood
# ------------------------------------------------------------

fit.truncated <- function(dat, interval, K, lambda) {
  
  Y <- dat$Y
  X <- dat$X
  Z <- dat$Z
  
  B <- make.spline.basis(Z, K)
  
  q <- ncol(B)
  
  # Initial values
  init <- c(
    beta.true,
    rep(0, q)
  )
  
  # ----------------------------------------------------------
  # Negative penalized truncated log-likelihood
  # ----------------------------------------------------------
  
  objective <- function(par) {
    
    beta <- par[1]
    
    theta <- par[-1]
    
    eta <- beta * X + as.vector(B %*% theta)
    
    # Regression density
    log.f <- dnorm(
      Y,
      mean = eta,
      sd = sigma,
      log = TRUE
    )
    
    # --------------------------------------------------------
    # Truncation probability
    #
    # P(a <= Y <= b | X,Z)
    # --------------------------------------------------------
    
    C <- pnorm(
      interval[2],
      mean = eta,
      sd = sigma
    ) -
      pnorm(
        interval[1],
        mean = eta,
        sd = sigma
      )
    
    # Numerical protection
    C <- pmax(C, 1e-10)
    
    # --------------------------------------------------------
    # Roughness penalty
    # --------------------------------------------------------
    
    penalty <- lambda * sum(diff(theta, differences = 2)^2)
    
    # Penalized negative log likelihood
    value <- -sum(log.f - log(C)) + penalty
    
    if (!is.finite(value)) {
      
      value <- 1e100
      
    }
    
    value
    
  }
  
  fit <- optim(
    par = init,
    fn = objective,
    method = "BFGS",
    control = list(
      maxit = 1000,
      reltol = 1e-8
    )
  )
  
  beta.hat <- fit$par[1]
  
  theta.hat <- fit$par[-1]
  
  list(
    beta = beta.hat,
    theta = theta.hat,
    convergence = fit$convergence,
    basis = B
  )
  
}

# ------------------------------------------------------------
# 8. Fit conventional untruncated likelihood
# ------------------------------------------------------------

fit.untruncated <- function(dat, K, lambda) {
  
  Y <- dat$Y
  X <- dat$X
  Z <- dat$Z
  
  B <- make.spline.basis(Z, K)
  
  q <- ncol(B)
  
  init <- c(
    beta.true,
    rep(0, q)
  )
  
  # ----------------------------------------------------------
  # Conventional untruncated likelihood
  # ----------------------------------------------------------
  
  objective <- function(par) {
    
    beta <- par[1]
    
    theta <- par[-1]
    
    eta <- beta * X + as.vector(B %*% theta)
    
    log.f <- dnorm(
      Y,
      mean = eta,
      sd = sigma,
      log = TRUE
    )
    
    penalty <- lambda *
      sum(diff(theta, differences = 2)^2)
    
    value <- -sum(log.f) + penalty
    
    if (!is.finite(value)) {
      
      value <- 1e100
      
    }
    
    value
    
  }
  
  fit <- optim(
    par = init,
    fn = objective,
    method = "BFGS",
    control = list(
      maxit = 1000,
      reltol = 1e-8
    )
  )
  
  beta.hat <- fit$par[1]
  
  theta.hat <- fit$par[-1]
  
  list(
    beta = beta.hat,
    theta = theta.hat,
    convergence = fit$convergence,
    basis = B
  )
  
}

# ------------------------------------------------------------
# 9. Prediction function
# ------------------------------------------------------------

predict.sieve <- function(fit, X, Z, K) {
  
  B <- make.spline.basis(Z, K)
  
  fit$beta * X +
    as.vector(B %*% fit$theta)
  
}

# ------------------------------------------------------------
# 10. L2 error of nonlinear component
# ------------------------------------------------------------

calculate.L2 <- function(fit, K) {
  
  z.grid <- seq(0, 1, length.out = 500)
  
  B.grid <- make.spline.basis(
    z.grid,
    K
  )
  
  m.hat <- as.vector(
    B.grid %*% fit$theta
  )
  
  m0 <- m.true(z.grid)
  
  sqrt(
    mean(
      (m.hat - m0)^2
    )
  )
  
}

# ------------------------------------------------------------
# 11. One Monte Carlo replication
# ------------------------------------------------------------

one.replication <- function(
    n,
    interval,
    regime
) {
  
  # Complete underlying sample
  complete.data <-
    generate.complete.data(n)
  
  # Apply truncation
  observed.data <-
    apply.truncation(
      complete.data,
      interval
    )
  
  # If too few observations remain,
  # repeat generation
  attempts <- 0
  
  while (
    nrow(observed.data) < max(50, 0.25 * n)
  ) {
    
    complete.data <-
      generate.complete.data(n)
    
    observed.data <-
      apply.truncation(
        complete.data,
        interval
      )
    
    attempts <- attempts + 1
    
    if (attempts > 100) {
      
      stop(
        "Too few observations after truncation."
      )
      
    }
    
  }
  
  # ----------------------------------------------------------
  # Sieve and penalty
  # ----------------------------------------------------------
  
  K <- get.K(n)
  
  lambda <-
    n^(-2 * r / (2 * r + 1))
  
  # ----------------------------------------------------------
  # Proposed estimator
  # ----------------------------------------------------------
  
  fit.T <- fit.truncated(
    observed.data,
    interval,
    K,
    lambda
  )
  
  # ----------------------------------------------------------
  # Conventional estimator
  # ----------------------------------------------------------
  
  fit.U <- fit.untruncated(
    observed.data,
    K,
    lambda
  )
  
  # ----------------------------------------------------------
  # Parametric estimation errors
  # ----------------------------------------------------------
  
  beta.error.T <-
    fit.T$beta - beta.true
  
  beta.error.U <-
    fit.U$beta - beta.true
  
  # ----------------------------------------------------------
  # Nonparametric L2 errors
  # ----------------------------------------------------------
  
  L2.T <-
    calculate.L2(
      fit.T,
      K
    )
  
  L2.U <-
    calculate.L2(
      fit.U,
      K
    )
  
  # ----------------------------------------------------------
  # Independent UNTRUNCATED test sample
  #
  # This evaluates recovery of the underlying population
  # regression relationship.
  # ----------------------------------------------------------
  
  test.data <-
    generate.complete.data(
      N.test
    )
  
  pred.T <-
    predict.sieve(
      fit.T,
      test.data$X,
      test.data$Z,
      K
    )
  
  pred.U <-
    predict.sieve(
      fit.U,
      test.data$X,
      test.data$Z,
      K
    )
  
  TestMSE.T <-
    mean(
      (test.data$Y - pred.T)^2
    )
  
  TestMSE.U <-
    mean(
      (test.data$Y - pred.U)^2
    )
  
  data.frame(
    
    n = n,
    
    Regime = regime,
    
    ObservedN = nrow(observed.data),
    
    Proposed_Beta_Error =
      beta.error.T,
    
    Untruncated_Beta_Error =
      beta.error.U,
    
    Proposed_L2 =
      L2.T,
    
    Untruncated_L2 =
      L2.U,
    
    Proposed_Test_MSE =
      TestMSE.T,
    
    Untruncated_Test_MSE =
      TestMSE.U
    
  )
  
}

# ------------------------------------------------------------
# 12. Monte Carlo simulation
# ------------------------------------------------------------

all.results <- list()

counter <- 1

for (n in n.values) {
  
  for (regime in names(truncation.regimes)) {
    
    interval <-
      truncation.regimes[[regime]]
    
    cat(
      "\nRunning:",
      "n =", n,
      "| Regime =", regime,
      "\n"
    )
    
    temp <- vector(
      "list",
      R
    )
    
    for (r.rep in 1:R) {
      
      temp[[r.rep]] <-
        one.replication(
          n = n,
          interval = interval,
          regime = regime
        )
      
    }
    
    all.results[[counter]] <-
      do.call(
        rbind,
        temp
      )
    
    counter <- counter + 1
    
  }
  
}

results <-
  do.call(
    rbind,
    all.results
  )

# ------------------------------------------------------------
# 13. Summarize Monte Carlo results
# ------------------------------------------------------------

summary.results <- do.call(
  rbind,
  lapply(
    split(
      results,
      list(
        results$n,
        results$Regime
      )
    ),
    function(d) {
      
      data.frame(
        
        n = unique(d$n),
        
        Regime =
          unique(d$Regime),
        
        Method = c(
          "Proposed Truncated",
          "Untruncated"
        ),
        
        Bias = c(
          mean(
            d$Proposed_Beta_Error
          ),
          mean(
            d$Untruncated_Beta_Error
          )
        ),
        
        RMSE_Beta = c(
          
          sqrt(
            mean(
              d$Proposed_Beta_Error^2
            )
          ),
          
          sqrt(
            mean(
              d$Untruncated_Beta_Error^2
            )
          )
          
        ),
        
        L2_Error = c(
          
          mean(
            d$Proposed_L2
          ),
          
          mean(
            d$Untruncated_L2
          )
          
        ),
        
        Test_MSE = c(
          
          mean(
            d$Proposed_Test_MSE
          ),
          
          mean(
            d$Untruncated_Test_MSE
          )
          
        )
        
      )
      
    }
  )
)

# ------------------------------------------------------------
# 14. Display final table
# ------------------------------------------------------------

print(
  summary.results,
  row.names = FALSE
)

# ------------------------------------------------------------
# 15. Save results
# ------------------------------------------------------------

write.csv(
  summary.results,
  "Truncated_vs_Untruncated_Comparison.csv",
  row.names = FALSE
)

