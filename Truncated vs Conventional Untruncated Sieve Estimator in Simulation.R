# ============================================================
# REAL DATA:
# Proposed Truncated vs Conventional Untruncated Estimator
# ============================================================

library(splines)

# ============================================================
# 1. DATA
# ============================================================

dat <- automotive_chip_reliability_dataset[, c(
  "Failure.Rate....",
  "Thermal.Stress",
  "Technology.Node",
  "Gate.Density",
  "Power.Consumption",
  "Mechanical.Stress",
  "Electrical.Density"
)]

# Remove missing observations
dat <- dat[complete.cases(dat), ]

# Rename variables
names(dat) <- c(
  "FailureRate",
  "ThermalStress",
  "TechnologyNode",
  "GateDensity",
  "PowerConsumption",
  "MechanicalStress",
  "ElectricalDensity"
)

n <- nrow(dat)

cat("Number of observations:", n, "\n")

cat("\nVariables used:\n")
print(names(dat))

cat("\nSummary statistics:\n")
print(summary(dat))


# ============================================================
# 2. VARIABLE DEFINITIONS
# ============================================================

Yvar <- "FailureRate"
Zvar <- "ThermalStress"

Xvars <- c(
  "TechnologyNode",
  "GateDensity",
  "PowerConsumption",
  "MechanicalStress",
  "ElectricalDensity"
)

Y <- dat[[Yvar]]
Z <- dat[[Zvar]]


# ============================================================
# 3. PARAMETRIC DESIGN MATRIX
# ============================================================
#
# The non-intercept covariates are standardized to avoid severe
# numerical conditioning problems.
#
# The final estimates are transformed back to the original
# measurement scale before reporting.
# ============================================================

X_raw <- dat[, Xvars, drop = FALSE]

X_scaled <- as.data.frame(
  scale(
    X_raw,
    center = TRUE,
    scale = TRUE
  )
)

X <- model.matrix(
  ~ .,
  data = X_scaled
)

p <- ncol(X)

cat("\nNumber of parametric coefficients:", p, "\n")
print(colnames(X))


# ============================================================
# 4. SPLINE BASIS
# ============================================================

K <- 5

boundary_knots <- range(
  Z,
  na.rm = TRUE
)

# Five internal knots
internal_knots <- quantile(
  Z,
  probs = seq(
    0,
    1,
    length.out = K + 2
  )[-c(1, K + 2)],
  names = FALSE
)

# IMPORTANT:
# intercept = FALSE avoids duplication of the intercept already
# present in X.
B <- bs(
  Z,
  knots = internal_knots,
  degree = 3,
  intercept = FALSE,
  Boundary.knots = boundary_knots
)

q <- ncol(B)

cat(
  "\nNumber of spline coefficients:",
  q,
  "\n"
)


# ============================================================
# 5. DESIGN MATRIX
# ============================================================

D <- cbind(
  X,
  B
)

k <- ncol(D)

cat(
  "Total number of coefficients:",
  k,
  "\n"
)

cat(
  "Condition number of design matrix:",
  kappa(D),
  "\n"
)


# ============================================================
# 6. PENALTY MATRIX
# ============================================================

penalty_matrix <- function(q) {
  
  P <- matrix(
    0,
    nrow = q,
    ncol = q
  )
  
  if (q >= 3) {
    
    for (j in 2:(q - 1)) {
      
      P[j, j] <-
        P[j, j] + 1
      
      P[j + 1, j + 1] <-
        P[j + 1, j + 1] + 1
      
      P[j, j + 1] <-
        P[j, j + 1] - 1
      
      P[j + 1, j] <-
        P[j + 1, j] - 1
    }
  }
  
  P
}

P <- penalty_matrix(q)


# ============================================================
# 7. NUMERICALLY STABLE SOLVER
# ============================================================

safe_solve <- function(A, b) {
  
  # Symmetrize
  A <- (A + t(A)) / 2
  
  # Small numerical ridge
  ridge <- 1e-8 * max(
    1,
    max(abs(diag(A)))
  )
  
  A_stable <- A +
    ridge * diag(nrow(A))
  
  solve(
    A_stable,
    b
  )
}


# ============================================================
# 8. UNTRUNCATED ESTIMATOR
# ============================================================

untruncated_fit <- function(
    y,
    x,
    b,
    lambda = 1
) {
  
  D <- cbind(
    x,
    b
  )
  
  p <- ncol(x)
  q <- ncol(b)
  k <- ncol(D)
  
  # Full penalty matrix
  Pfull <- matrix(
    0,
    nrow = k,
    ncol = k
  )
  
  Pfull[
    (p + 1):(p + q),
    (p + 1):(p + q)
  ] <- P
  
  A <- crossprod(D) +
    lambda * Pfull
  
  rhs <- crossprod(
    D,
    y
  )
  
  beta_hat <- safe_solve(
    A,
    rhs
  )
  
  fitted <- as.vector(
    D %*% beta_hat
  )
  
  residuals <- y - fitted
  
  df_res <- max(
    1,
    length(y) - k
  )
  
  sigma2 <- sum(
    residuals^2
  ) / df_res
  
  # Sandwich-type covariance approximation
  Ainv <- safe_solve(
    A,
    diag(k)
  )
  
  XtX <- crossprod(D)
  
  vcov_beta <- sigma2 *
    Ainv %*%
    XtX %*%
    Ainv
  
  vcov_beta <- (
    vcov_beta +
      t(vcov_beta)
  ) / 2
  
  list(
    coef = as.numeric(beta_hat),
    fitted = fitted,
    residuals = residuals,
    sigma = sqrt(sigma2),
    vcov = vcov_beta,
    convergence = 0
  )
}


# ============================================================
# 9. FIT CONVENTIONAL UNTRUNCATED MODEL
# ============================================================

fit_U <- untruncated_fit(
  y = Y,
  x = X,
  b = B,
  lambda = 1
)

cat(
  "\nUntruncated estimator successfully fitted.\n"
)

cat(
  "Estimated sigma:",
  fit_U$sigma,
  "\n"
)


# ============================================================
# 10. TRUNCATION LIMITS
# ============================================================
#
# IMPORTANT:
# DO NOT use min(Y) and max(Y) here.
#
# Insert the truncation limits specified in the manuscript.
# ============================================================

# ------------------------------------------------------------
# Truncation limits
# ------------------------------------------------------------
# Based on the observed FailureRate range in the dataset

aY <- 1.01113
bY <- 9.959324

cat("\nSpecified truncation range:\n")
cat("[", aY, ", ", bY, "]\n", sep = "")

cat(
  "\nSpecified truncation range:\n"
)

cat(
  "[",
  aY,
  ", ",
  bY,
  "]\n",
  sep = ""
)


# ============================================================
# 11. STABLE LOG DIFFERENCE OF NORMAL CDFs
# ============================================================

log_pnorm_difference <- function(
    upper,
    lower
) {
  
  # Direct calculation first
  C <- pnorm(
    upper
  ) -
    pnorm(
      lower
    )
  
  # Prevent log(0)
  C <- pmax(
    C,
    1e-12
  )
  
  log(C)
}


# ============================================================
# 12. PROPOSED TRUNCATED ESTIMATOR
# ============================================================

truncated_fit <- function(
    y,
    x,
    b,
    aY,
    bY,
    lambda = 1
) {
  
  D <- cbind(
    x,
    b
  )
  
  p <- ncol(x)
  q <- ncol(b)
  k <- ncol(D)
  
  # Full penalty matrix
  Pfull <- matrix(
    0,
    nrow = k,
    ncol = k
  )
  
  Pfull[
    (p + 1):(p + q),
    (p + 1):(p + q)
  ] <- P
  
  # ----------------------------------------------------------
  # Initial estimator
  # ----------------------------------------------------------
  
  A <- crossprod(D) +
    lambda * Pfull
  
  rhs <- crossprod(
    D,
    y
  )
  
  theta0 <- safe_solve(
    A,
    rhs
  )
  
  mu0 <- as.vector(
    D %*% theta0
  )
  
  sigma0 <- sqrt(
    mean(
      (y - mu0)^2
    )
  )
  
  sigma0 <- max(
    sigma0,
    1e-4
  )
  
  # ----------------------------------------------------------
  # Parameter vector
  # theta = coefficients
  # log_sigma = log(sigma)
  # ----------------------------------------------------------
  
  par0 <- c(
    theta0,
    log(sigma0)
  )
  
  # ----------------------------------------------------------
  # Objective function
  # ----------------------------------------------------------
  
  objective <- function(par) {
    
    theta <- par[
      1:k
    ]
    
    log_sigma <- par[
      k + 1
    ]
    
    sigma <- exp(
      log_sigma
    )
    
    # Avoid numerical overflow
    sigma <- max(
      sigma,
      1e-6
    )
    
    mu <- as.vector(
      D %*% theta
    )
    
    upper <- (
      bY - mu
    ) / sigma
    
    lower <- (
      aY - mu
    ) / sigma
    
    # Truncation probability
    C <- pnorm(
      upper
    ) -
      pnorm(
        lower
      )
    
    # Numerical protection
    C <- pmax(
      C,
      1e-12
    )
    
    loglik <- sum(
      dnorm(
        y,
        mean = mu,
        sd = sigma,
        log = TRUE
      ) -
        log(C)
    )
    
    # Spline penalty
    spline_coef <- theta[
      (p + 1):(p + q)
    ]
    
    penalty <- lambda *
      as.numeric(
        crossprod(
          spline_coef,
          P %*% spline_coef
        )
      )
    
    -(loglik - penalty)
  }
  
  
  # ----------------------------------------------------------
  # Optimization
  # ----------------------------------------------------------
  
  opt <- optim(
    par = par0,
    fn = objective,
    method = "BFGS",
    hessian = TRUE,
    control = list(
      maxit = 3000,
      reltol = 1e-8
    )
  )
  
  
  # ----------------------------------------------------------
  # Check optimization
  # ----------------------------------------------------------
  
  if (
    opt$convergence != 0 ||
    !is.finite(opt$value)
  ) {
    
    warning(
      "Truncated optimization did not fully converge."
    )
  }
  
  
  # ----------------------------------------------------------
  # Extract estimates
  # ----------------------------------------------------------
  
  theta_hat <- opt$par[
    1:k
  ]
  
  sigma_hat <- exp(
    opt$par[
      k + 1
    ]
  )
  
  fitted <- as.vector(
    D %*% theta_hat
  )
  
  # ----------------------------------------------------------
  # Covariance matrix
  # ----------------------------------------------------------
  
  vcov_full <- tryCatch(
    
    solve(
      (
        opt$hessian +
          t(opt$hessian)
      ) / 2
    ),
    
    error = function(e) NULL
  )
  
  if (
    !is.null(vcov_full)
  ) {
    
    vcov_full <- (
      vcov_full +
        t(vcov_full)
    ) / 2
  }
  
  list(
    coef = as.numeric(theta_hat),
    sigma = sigma_hat,
    fitted = fitted,
    convergence = opt$convergence,
    objective = opt$value,
    hessian = opt$hessian,
    vcov = vcov_full
  )
}


# ============================================================
# 13. FIT PROPOSED MODEL
# ============================================================

fit_T <- truncated_fit(
  y = Y,
  x = X,
  b = B,
  aY = aY,
  bY = bY,
  lambda = 1
)

cat(
  "\nProposed estimator convergence code:",
  fit_T$convergence,
  "\n"
)

cat(
  "Proposed estimator sigma:",
  fit_T$sigma,
  "\n"
)


# ============================================================
# 14. PARAMETRIC COEFFICIENT NAMES
# ============================================================

coef_names <- c(
  "Intercept",
  "TechnologyNode",
  "GateDensity",
  "PowerConsumption",
  "MechanicalStress",
  "ElectricalDensity"
)


# ============================================================
# 15. CONVERT STANDARDIZED COEFFICIENTS BACK TO ORIGINAL SCALE
# ============================================================

x_center <- attr(
  scale(X_raw),
  "scaled:center"
)

x_scale <- attr(
  scale(X_raw),
  "scaled:scale"
)

# The values above need to be calculated explicitly
# to avoid repeated scale() calls.

x_center <- sapply(
  X_raw,
  mean
)

x_scale <- sapply(
  X_raw,
  sd
)


# ------------------------------------------------------------
# Transformation function
# ------------------------------------------------------------

back_transform <- function(
    beta_scaled,
    x_center,
    x_scale
) {
  
  beta_original <- numeric(
    length(beta_scaled)
  )
  
  names(beta_original) <- coef_names
  
  # Slopes
  beta_original[
    -1
  ] <- beta_scaled[
    -1
  ] / x_scale
  
  # Intercept
  beta_original[
    1
  ] <- beta_scaled[
    1
  ] -
    sum(
      beta_scaled[-1] *
        x_center /
        x_scale
    )
  
  beta_original
}


# ============================================================
# 16. PARAMETRIC ESTIMATES
# ============================================================

beta_T_scaled <- fit_T$coef[
  1:p
]

beta_U_scaled <- fit_U$coef[
  1:p
]


beta_T <- back_transform(
  beta_T_scaled,
  x_center,
  x_scale
)

beta_U <- back_transform(
  beta_U_scaled,
  x_center,
  x_scale
)


# ============================================================
# 17. STANDARD ERRORS
# ============================================================

se_T_scaled <- rep(
  NA_real_,
  p
)

if (
  !is.null(fit_T$vcov)
) {
  
  d <- diag(
    fit_T$vcov
  )
  
  d <- pmax(
    d,
    0
  )
  
  se_T_scaled <- sqrt(
    d[1:p]
  )
}


se_U_scaled <- rep(
  NA_real_,
  p
)

if (
  !is.null(fit_U$vcov)
) {
  
  d <- diag(
    fit_U$vcov
  )
  
  d <- pmax(
    d,
    0
  )
  
  se_U_scaled <- sqrt(
    d[1:p]
  )
}


# Transform slope SEs back to original scale

se_T <- c(
  se_T_scaled[1],
  se_T_scaled[-1] / x_scale
)

se_U <- c(
  se_U_scaled[1],
  se_U_scaled[-1] / x_scale
)


# ============================================================
# 18. COEFFICIENT COMPARISON TABLE
# ============================================================

coef_table <- data.frame(
  
  Variable = coef_names,
  
  Proposed_Estimate =
    as.numeric(beta_T),
  
  Proposed_SE =
    as.numeric(se_T),
  
  Untruncated_Estimate =
    as.numeric(beta_U),
  
  Untruncated_SE =
    as.numeric(se_U),
  
  row.names = NULL
)


cat(
  "\n============================================================\n"
)

cat(
  "PARAMETRIC COEFFICIENT COMPARISON\n"
)

cat(
  "============================================================\n"
)

print(
  coef_table,
  digits = 6,
  row.names = FALSE
)


# ============================================================
# 19. COMMON SPLINE GRID
# ============================================================

z_grid <- seq(
  min(Z),
  max(Z),
  length.out = 200
)

B_grid <- bs(
  z_grid,
  knots = internal_knots,
  degree = 3,
  intercept = FALSE,
  Boundary.knots = boundary_knots
)


# ============================================================
# 20. SMOOTH COMPONENT
# ============================================================

theta_T_smooth <- fit_T$coef[
  (p + 1):(p + q)
]

theta_U_smooth <- fit_U$coef[
  (p + 1):(p + q)
]


m_T <- as.vector(
  B_grid %*%
    theta_T_smooth
)

m_U <- as.vector(
  B_grid %*%
    theta_U_smooth
)


smooth_comparison <- data.frame(
  
  ThermalStress = z_grid,
  
  Proposed_Truncated = m_T,
  
  Untruncated = m_U
)


cat(
  "\n============================================================\n"
)

cat(
  "ESTIMATED NONLINEAR FUNCTION\n"
)

cat(
  "============================================================\n"
)

print(
  head(
    smooth_comparison
  ),
  row.names = FALSE
)


# ============================================================
# 21. PLOT
# ============================================================

if (
  all(
    is.finite(m_T)
  ) &&
  all(
    is.finite(m_U)
  )
) {
  
  pdf(
    "Real_Data_Smooth_Comparison.pdf",
    width = 8,
    height = 5.5
  )
  
  plot(
    z_grid,
    m_T,
    type = "l",
    lwd = 2,
    xlab = "Thermal Stress",
    ylab = "Estimated smooth effect",
    main = "Estimated Nonlinear Thermal Stress Effect"
  )
  
  lines(
    z_grid,
    m_U,
    lwd = 2,
    lty = 2
  )
  
  rug(
    Z
  )
  
  legend(
    "topright",
    legend = c(
      "Proposed truncated",
      "Conventional untruncated"
    ),
    lty = c(
      1,
      2
    ),
    lwd = 2,
    bty = "n"
  )
  
  dev.off()
  
} else {
  
  warning(
    "Smooth estimates contain non-finite values; plot was not generated."
  )
}


# ============================================================
# 22. 5-FOLD CROSS-VALIDATION
# ============================================================

set.seed(123)

Kfold <- 5

fold_id <- sample(
  rep(
    1:Kfold,
    length.out = n
  )
)

cv_results <- data.frame()


for (
  fold in 1:Kfold
) {
  
  cat(
    "\nRunning fold:",
    fold,
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Training and test data
  # ----------------------------------------------------------
  
  train <- dat[
    fold_id != fold,
    ,
    drop = FALSE
  ]
  
  test <- dat[
    fold_id == fold,
    ,
    drop = FALSE
  ]
  
  
  # ----------------------------------------------------------
  # Standardization using TRAINING data only
  # ----------------------------------------------------------
  
  train_raw <- train[
    Xvars
  ]
  
  test_raw <- test[
    Xvars
  ]
  
  train_center <- sapply(
    train_raw,
    mean
  )
  
  train_scale <- sapply(
    train_raw,
    sd
  )
  
  # Protect against zero variance
  train_scale[
    train_scale == 0
  ] <- 1
  
  
  train_scaled <- as.data.frame(
    scale(
      train_raw,
      center = train_center,
      scale = train_scale
    )
  )
  
  test_scaled <- as.data.frame(
    scale(
      test_raw,
      center = train_center,
      scale = train_scale
    )
  )
  
  
  # ----------------------------------------------------------
  # Parametric matrices
  # ----------------------------------------------------------
  
  X_train <- model.matrix(
    ~ .,
    data = train_scaled
  )
  
  X_test <- model.matrix(
    ~ .,
    data = test_scaled
  )
  
  
  # ----------------------------------------------------------
  # Common spline basis
  # ----------------------------------------------------------
  
  B_train <- bs(
    train[[Zvar]],
    knots = internal_knots,
    degree = 3,
    intercept = FALSE,
    Boundary.knots = boundary_knots
  )
  
  B_test <- bs(
    test[[Zvar]],
    knots = internal_knots,
    degree = 3,
    intercept = FALSE,
    Boundary.knots = boundary_knots
  )
  
  
  # ----------------------------------------------------------
  # Responses
  # ----------------------------------------------------------
  
  y_train <- train[[Yvar]]
  y_test <- test[[Yvar]]
  
  
  # ----------------------------------------------------------
  # Fit proposed estimator
  # ----------------------------------------------------------
  
  fit_T_cv <- tryCatch(
    
    truncated_fit(
      y = y_train,
      x = X_train,
      b = B_train,
      aY = aY,
      bY = bY,
      lambda = 1
    ),
    
    error = function(e) {
      
      warning(
        paste(
          "Truncated estimator failed in fold",
          fold,
          ":",
          e$message
        )
      )
      
      NULL
    }
  )
  
  
  # ----------------------------------------------------------
  # Fit conventional estimator
  # ----------------------------------------------------------
  
  fit_U_cv <- tryCatch(
    
    untruncated_fit(
      y = y_train,
      x = X_train,
      b = B_train,
      lambda = 1
    ),
    
    error = function(e) {
      
      warning(
        paste(
          "Untruncated estimator failed in fold",
          fold,
          ":",
          e$message
        )
      )
      
      NULL
    }
  )
  
  
  # ----------------------------------------------------------
  # Check fitting
  # ----------------------------------------------------------
  
  if (
    is.null(fit_T_cv) ||
    is.null(fit_U_cv)
  ) {
    
    next
  }
  
  
  # ----------------------------------------------------------
  # Predictions
  # ----------------------------------------------------------
  
  D_test <- cbind(
    X_test,
    B_test
  )
  
  
  pred_T <- as.vector(
    D_test %*%
      fit_T_cv$coef
  )
  
  pred_U <- as.vector(
    D_test %*%
      fit_U_cv$coef
  )
  
  
  # ----------------------------------------------------------
  # Check predictions
  # ----------------------------------------------------------
  
  if (
    any(
      !is.finite(pred_T)
    ) ||
    any(
      !is.finite(pred_U)
    )
  ) {
    
    warning(
      paste(
        "Non-finite prediction in fold",
        fold
      )
    )
    
    next
  }
  
  
  # ----------------------------------------------------------
  # RMSE
  # ----------------------------------------------------------
  
  RMSE_T <- sqrt(
    mean(
      (
        y_test -
          pred_T
      )^2
    )
  )
  
  RMSE_U <- sqrt(
    mean(
      (
        y_test -
          pred_U
      )^2
    )
  )
  
  
  # ----------------------------------------------------------
  # MAE
  # ----------------------------------------------------------
  
  MAE_T <- mean(
    abs(
      y_test -
        pred_T
    )
  )
  
  MAE_U <- mean(
    abs(
      y_test -
        pred_U
    )
  )
  
  
  # ----------------------------------------------------------
  # R-squared
  # ----------------------------------------------------------
  
  SST <- sum(
    (
      y_test -
        mean(y_train)
    )^2
  )
  
  R2_T <- 1 -
    sum(
      (
        y_test -
          pred_T
      )^2
    ) /
    SST
  
  R2_U <- 1 -
    sum(
      (
        y_test -
          pred_U
      )^2
    ) /
    SST
  
  
  # ----------------------------------------------------------
  # Store results
  # ----------------------------------------------------------
  
  cv_results <- rbind(
    
    cv_results,
    
    data.frame(
      
      Fold = fold,
      
      Method = c(
        "Proposed Truncated",
        "Untruncated"
      ),
      
      RMSE = c(
        RMSE_T,
        RMSE_U
      ),
      
      MAE = c(
        MAE_T,
        MAE_U
      ),
      
      R2 = c(
        R2_T,
        R2_U
      )
    )
  )
}


# ============================================================
# 23. CROSS-VALIDATION RESULTS
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "5-FOLD CROSS-VALIDATION RESULTS\n"
)

cat(
  "============================================================\n"
)

print(
  cv_results,
  row.names = FALSE
)


# ============================================================
# 24. CROSS-VALIDATED SUMMARY
# ============================================================

if (
  nrow(cv_results) > 0
) {
  
  cv_summary <- do.call(
    rbind,
    lapply(
      split(
        cv_results,
        cv_results$Method
      ),
      function(d) {
        
        data.frame(
          
          Method =
            unique(
              d$Method
            ),
          
          CV_RMSE =
            mean(
              d$RMSE,
              na.rm = TRUE
            ),
          
          SD_RMSE =
            sd(
              d$RMSE,
              na.rm = TRUE
            ),
          
          CV_MAE =
            mean(
              d$MAE,
              na.rm = TRUE
            ),
          
          SD_MAE =
            sd(
              d$MAE,
              na.rm = TRUE
            ),
          
          CV_R2 =
            mean(
              d$R2,
              na.rm = TRUE
            ),
          
          SD_R2 =
            sd(
              d$R2,
              na.rm = TRUE
            )
        )
      }
    )
  )
  
  
  cat(
    "\n============================================================\n"
  )
  
  cat(
    "CROSS-VALIDATED PREDICTIVE PERFORMANCE\n"
  )
  
  cat(
    "============================================================\n"
  )
  
  print(
    cv_summary,
    digits = 6,
    row.names = FALSE
  )
  
  
  # ==========================================================
  # 25. PREDICTIVE IMPROVEMENT
  # ==========================================================
  
  rmse_T <- cv_summary$CV_RMSE[
    cv_summary$Method ==
      "Proposed Truncated"
  ]
  
  rmse_U <- cv_summary$CV_RMSE[
    cv_summary$Method ==
      "Untruncated"
  ]
  
  mae_T <- cv_summary$CV_MAE[
    cv_summary$Method ==
      "Proposed Truncated"
  ]
  
  mae_U <- cv_summary$CV_MAE[
    cv_summary$Method ==
      "Untruncated"
  ]
  
  
  if (
    length(rmse_T) == 1 &&
    length(rmse_U) == 1 &&
    length(mae_T) == 1 &&
    length(mae_U) == 1
  ) {
    
    rmse_improvement <- 100 *
      (
        rmse_U -
          rmse_T
      ) /
      rmse_U
    
    mae_improvement <- 100 *
      (
        mae_U -
          mae_T
      ) /
      mae_U
    
    
    cat(
      "\n============================================================\n"
    )
    
    cat(
      "PREDICTIVE IMPROVEMENT\n"
    )
    
    cat(
      "============================================================\n"
    )
    
    cat(
      "RMSE improvement (%):",
      round(
        rmse_improvement,
        2
      ),
      "\n"
    )
    
    cat(
      "MAE improvement (%):",
      round(
        mae_improvement,
        2
      ),
      "\n"
    )
    
  }
  
} else {
  
  warning(
    "No valid cross-validation results were produced."
  )
}


#################################################################################


# ============================================================
# GENERATE SMOOTH FUNCTION COMPARISON PLOT
# ============================================================

# Grid of Thermal Stress values
z_grid <- seq(
  min(Z, na.rm = TRUE),
  max(Z, na.rm = TRUE),
  length.out = 200
)

# ------------------------------------------------------------
# Construct spline basis on the same basis used for estimation
# ------------------------------------------------------------

B_grid <- bs(
  z_grid,
  knots = internal_knots,
  degree = 3,
  intercept = FALSE,
  Boundary.knots = boundary_knots
)

# ------------------------------------------------------------
# Extract spline coefficients
# ------------------------------------------------------------

theta_T_smooth <- fit_T$coef[
  (p + 1):(p + q)
]

theta_U_smooth <- fit_U$coef[
  (p + 1):(p + q)
]

# ------------------------------------------------------------
# Estimated smooth functions
# ------------------------------------------------------------

m_T <- as.vector(
  B_grid %*% theta_T_smooth
)

m_U <- as.vector(
  B_grid %*% theta_U_smooth
)

# ------------------------------------------------------------
# Check dimensions and finite values
# ------------------------------------------------------------

cat("Length of z_grid:", length(z_grid), "\n")
cat("Length of m_T:", length(m_T), "\n")
cat("Length of m_U:", length(m_U), "\n")

cat(
  "Proposed estimator - finite values:",
  all(is.finite(m_T)),
  "\n"
)

cat(
  "Untruncated estimator - finite values:",
  all(is.finite(m_U)),
  "\n"
)

# ------------------------------------------------------------
# Generate PDF
# ------------------------------------------------------------

if (
  all(is.finite(z_grid)) &&
  all(is.finite(m_T)) &&
  all(is.finite(m_U))
) {
  
  pdf(
    file = "Real_Data_Smooth_Comparison.pdf",
    width = 8,
    height = 5.5
  )
  
  plot(
    z_grid,
    m_T,
    type = "l",
    lwd = 2,
    xlab = "Thermal Stress",
    ylab = "Estimated smooth effect",
    main = "Estimated Nonlinear Thermal Stress Effect"
  )
  
  lines(
    z_grid,
    m_U,
    lwd = 2,
    lty = 2
  )
  
  rug(
    Z,
    ticksize = 0.03
  )
  
  legend(
    "topleft",
    legend = c(
      "Proposed truncated",
      "Conventional untruncated"
    ),
    lty = c(1, 2),
    lwd = 2,
    bty = "n"
  )
  
  dev.off()
  
  cat(
    "\nPDF successfully generated:\n",
    normalizePath(
      "Real_Data_Smooth_Comparison.pdf"
    ),
    "\n"
  )
  
} else {
  
  stop(
    "The smooth estimates contain non-finite values. ",
    "PDF was not generated."
  )
  
}

