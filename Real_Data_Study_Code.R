# ============================================================
# Real Data Application: Automotive Chip Reliability Dataset
# ============================================================

# Install required package (run only once if not installed)
install.packages("ggplot2", repos = "https://cloud.r-project.org")

# Load package
library(ggplot2)

# ------------------------------------------------------------
# 1. Read the dataset
# ------------------------------------------------------------

data <- automotive_chip_reliability_dataset
names(data)

# ============================================================
# Select variables used in the proposed model
# ============================================================

dat <- data[, c(
  "Failure.Rate....",
  "Thermal.Stress",
  "Technology.Node",
  "Gate.Density",
  "Power.Consumption",
  "Mechanical.Stress",
  "Electrical.Density"
)]

# Remove missing observations
dat <- na.omit(dat)


# ============================================================
# Rename variables for convenient analysis
# ============================================================

names(dat) <- c(
  "FailureRate",
  "ThermalStress",
  "TechnologyNode",
  "GateDensity",
  "PowerConsumption",
  "MechanicalStress",
  "ElectricalDensity"
)


# ============================================================
# Basic inspection
# ============================================================

dim(dat)
names(dat)
str(dat)
summary(dat)

# Number of observations
nrow(dat)


# ------------------------------------------------------------
# 3. Rename variables for convenient analysis
# ------------------------------------------------------------

names(dat) <- c(
  "FailureRate",
  "ThermalStress",
  "TechnologyNode",
  "GateDensity",
  "PowerConsumption",
  "MechanicalStress",
  "ElectricalDensity"
)

str(dat)
summary(dat)


# ============================================================
# 4. DATA VISUALIZATION
# ============================================================


# ------------------------------------------------------------
# 4.1 Distribution of Failure Rate
# ------------------------------------------------------------

p1 <- ggplot(dat, aes(x = FailureRate)) +
  geom_histogram(bins = 15, 
                 fill = "steelblue",
                 colour = "black") +
  labs(
    title = "Distribution of Failure Rate",
    x = "Failure Rate (%)",
    y = "Frequency"
  ) +
  theme_minimal()

print(p1)


# ------------------------------------------------------------
# 4.2 Density plot of Failure Rate
# ------------------------------------------------------------

p2 <- ggplot(dat, aes(x = FailureRate)) +
  geom_density(fill = "steelblue",
               alpha = 0.4,
               colour = "black") +
  labs(
    title = "Density of Failure Rate",
    x = "Failure Rate (%)",
    y = "Density"
  ) +
  theme_minimal()

print(p2)


# ------------------------------------------------------------
# 4.3 Failure Rate versus Thermal Stress
# ------------------------------------------------------------

p3 <- ggplot(dat, aes(x = ThermalStress,
                      y = FailureRate)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm",
              se = TRUE) +
  labs(
    title = "Failure Rate versus Thermal Stress",
    x = "Thermal Stress",
    y = "Failure Rate (%)"
  ) +
  theme_minimal()

print(p3)


# ------------------------------------------------------------
# 4.4 Nonlinear relationship: LOESS curve
# ------------------------------------------------------------

p4 <- ggplot(dat, aes(x = ThermalStress,
                      y = FailureRate)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess",
              formula = y ~ x,
              se = TRUE) +
  labs(
    title = "Nonlinear Relationship between Failure Rate and Thermal Stress",
    x = "Thermal Stress",
    y = "Failure Rate (%)"
  ) +
  theme_minimal()

print(p4)


# ------------------------------------------------------------
# 4.5 Linear and nonlinear curves together
# ------------------------------------------------------------

p5 <- ggplot(dat, aes(x = ThermalStress,
                      y = FailureRate)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm",
              se = FALSE,
              linetype = "dashed") +
  geom_smooth(method = "loess",
              se = FALSE) +
  labs(
    title = "Linear and Nonlinear Fits",
    x = "Thermal Stress",
    y = "Failure Rate (%)"
  ) +
  theme_minimal()

print(p5)


# ------------------------------------------------------------
# 4.6 Failure Rate versus Mechanical Stress
# ------------------------------------------------------------

p6 <- ggplot(dat, aes(x = MechanicalStress,
                      y = FailureRate)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess",
              se = TRUE) +
  labs(
    title = "Failure Rate versus Mechanical Stress",
    x = "Mechanical Stress",
    y = "Failure Rate (%)"
  ) +
  theme_minimal()

print(p6)


# ------------------------------------------------------------
# 4.7 Failure Rate versus Power Consumption
# ------------------------------------------------------------

p7 <- ggplot(dat, aes(x = PowerConsumption,
                      y = FailureRate)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess",
              se = TRUE) +
  labs(
    title = "Failure Rate versus Power Consumption",
    x = "Power Consumption",
    y = "Failure Rate (%)"
  ) +
  theme_minimal()

print(p7)


# ------------------------------------------------------------
# 4.8 Failure Rate versus Gate Density
# ------------------------------------------------------------

p8 <- ggplot(dat, aes(x = GateDensity,
                      y = FailureRate)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess",
              se = TRUE) +
  labs(
    title = "Failure Rate versus Gate Density",
    x = "Gate Density",
    y = "Failure Rate (%)"
  ) +
  theme_minimal()

print(p8)


# ------------------------------------------------------------
# 4.9 Failure Rate versus Electrical Density
# ------------------------------------------------------------

p9 <- ggplot(dat, aes(x = ElectricalDensity,
                      y = FailureRate)) +
  geom_point(size = 2) +
  geom_smooth(method = "loess",
              se = TRUE) +
  labs(
    title = "Failure Rate versus Electrical Density",
    x = "Electrical Density",
    y = "Failure Rate (%)"
  ) +
  theme_minimal()

print(p9)


# ============================================================
# 5. CORRELATION ANALYSIS
# ============================================================

cor_data <- dat[, c(
  "FailureRate",
  "ThermalStress",
  "TechnologyNode",
  "GateDensity",
  "PowerConsumption",
  "MechanicalStress",
  "ElectricalDensity"
)]

# Check the selected data
str(cor_data)

# Correlation matrix
cor_matrix <- cor(
  cor_data,
  use = "complete.obs",
  method = "pearson"
)

# Display correlation matrix
round(cor_matrix, 3)


# ------------------------------------------------------------
# 5.1 Correlation plot
# ------------------------------------------------------------

cor_matrix <- cor(
  cor_data,
  use = "complete.obs",
  method = "pearson"
)

# Display correlation matrix
round(cor_matrix, 2)

# Save correlation plot
pdf("Correlation_Plot.pdf", width = 9, height = 8)

par(mar = c(10, 10, 4, 2))

image(
  1:ncol(cor_matrix),
  1:nrow(cor_matrix),
  t(cor_matrix[nrow(cor_matrix):1, ]),
  axes = FALSE,
  col = gray.colors(100),
  main = "Correlation Matrix"
)

axis(
  1,
  at = 1:ncol(cor_matrix),
  labels = colnames(cor_matrix),
  las = 2
)

axis(
  2,
  at = 1:nrow(cor_matrix),
  labels = rev(rownames(cor_matrix)),
  las = 2
)

# Add correlation coefficients
for (i in 1:nrow(cor_matrix)) {
  for (j in 1:ncol(cor_matrix)) {
    text(
      j,
      nrow(cor_matrix) - i + 1,
      labels = round(cor_matrix[i, j], 2),
      cex = 0.8
    )
  }
}

dev.off()

# ============================================================
# 6. MULTIVARIATE VISUALIZATION
# ============================================================

pairs(
  cor_data,
  pch = 19,
  main = "Pairwise Relationships among Model Variables"
)


# ============================================================
# 7. BASIC DESCRIPTIVE MEASURES
# ============================================================

descriptive <- data.frame(
  Variable = names(cor_data),
  Mean = sapply(cor_data, mean, na.rm = TRUE),
  SD = sapply(cor_data, sd, na.rm = TRUE),
  Minimum = sapply(cor_data, min, na.rm = TRUE),
  Median = sapply(cor_data, median, na.rm = TRUE),
  Maximum = sapply(cor_data, max, na.rm = TRUE)
)

descriptive


# ============================================================
# Real Data: Partially Linear Spline Model
# ============================================================

# Required package
library(mgcv)

# Fit partially linear model
# ThermalStress = nonlinear component
# Other variables = linear parametric component

fit_pl <- gam(
  FailureRate ~
    TechnologyNode +
    GateDensity +
    PowerConsumption +
    MechanicalStress +
    ElectricalDensity +
    s(ThermalStress, k = 5, bs = "tp"),
  data = dat,
  method = "REML"
)

summary(fit_pl)
# Parametric coefficient estimates
coef_table <- summary(fit_pl)$p.table

coef_table
# Smooth-term summary
smooth_table <- summary(fit_pl)$s.table

smooth_table
pdf("Estimated_Nonlinear_Thermal_Stress_Effect.pdf",
    width = 7, height = 5)

plot(
  fit_pl,
  select = 1,
  shade = TRUE,
  residuals = TRUE,
  rug = TRUE,
  xlab = "Thermal Stress",
  ylab = "Estimated smooth effect",
  main = "Estimated Nonlinear Effect of Thermal Stress"
)

dev.off()
# Fitted values
fitted_values <- fitted(fit_pl)

# Residuals
residuals_pl <- residuals(fit_pl)

# RMSE
RMSE <- sqrt(mean(residuals_pl^2, na.rm = TRUE))

# MAE
MAE <- mean(abs(residuals_pl), na.rm = TRUE)

# R-squared
SSE <- sum(residuals_pl^2, na.rm = TRUE)
SST <- sum(
  (dat$FailureRate - mean(dat$FailureRate, na.rm = TRUE))^2,
  na.rm = TRUE
)

R2 <- 1 - SSE / SST

performance <- data.frame(
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2
)

performance
pdf("Residual_Diagnostics.pdf",
    width = 8, height = 6)

par(mfrow = c(2, 2))

plot(
  fitted_values,
  residuals_pl,
  pch = 19,
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs Fitted"
)

abline(h = 0, lty = 2)

plot(
  dat$ThermalStress,
  residuals_pl,
  pch = 19,
  xlab = "Thermal Stress",
  ylab = "Residuals",
  main = "Residuals vs Thermal Stress"
)

abline(h = 0, lty = 2)

hist(
  residuals_pl,
  breaks = 12,
  main = "Distribution of Residuals",
  xlab = "Residuals"
)

qqnorm(
  residuals_pl,
  main = "Normal Q-Q Plot"
)

qqline(residuals_pl)

dev.off()
# Fully linear model
fit_linear <- lm(
  FailureRate ~
    TechnologyNode +
    GateDensity +
    PowerConsumption +
    MechanicalStress +
    ElectricalDensity +
    ThermalStress,
  data = dat
)

# Compare AIC
AIC(fit_linear, fit_pl)
pred_linear <- fitted(fit_linear)
res_linear <- residuals(fit_linear)

RMSE_linear <- sqrt(mean(res_linear^2))
MAE_linear <- mean(abs(res_linear))

SSE_linear <- sum(res_linear^2)
R2_linear <- 1 -
  SSE_linear /
  sum((dat$FailureRate -
         mean(dat$FailureRate))^2)

linear_performance <- data.frame(
  Model = "Fully Linear",
  RMSE = RMSE_linear,
  MAE = MAE_linear,
  R2 = R2_linear
)

pl_performance <- data.frame(
  Model = "Partially Linear Spline",
  RMSE = RMSE,
  MAE = MAE,
  R2 = R2
)

comparison <- rbind(
  linear_performance,
  pl_performance
)

comparison

################################################################################

# ============================================================
# Estimated smooth effect of Thermal Stress
# ============================================================

library(mgcv)

# Fit the partially nonlinear model
fit_pl <- gam(
  FailureRate ~
    TechnologyNode +
    GateDensity +
    PowerConsumption +
    MechanicalStress +
    ElectricalDensity +
    s(ThermalStress, k = 5, bs = "tp"),
  data = dat,
  method = "REML"
)

# Save the estimated smooth effect
pdf("Estimated_Nonlinear_Thermal_Stress_Effect.pdf",
    width = 7, height = 5)

plot(
  fit_pl,
  select = 1,
  shade = TRUE,
  residuals = TRUE,
  rug = TRUE,
  xlab = "Thermal Stress",
  ylab = "Estimated smooth effect",
  main = "Estimated Effect of Thermal Stress"
)

dev.off()

################################################################################

# ============================================================
# Residual Diagnostics for the Partially Linear Spline Model
# ============================================================

# Fitted values and residuals
fitted_values <- fitted(fit_pl)
residuals_pl <- residuals(fit_pl)

# Save residual diagnostic plots
pdf("Residual_Diagnostics.pdf",
    width = 8, height = 6)

par(mfrow = c(2, 2))

# 1. Residuals versus fitted values
plot(
  fitted_values,
  residuals_pl,
  pch = 19,
  xlab = "Fitted values",
  ylab = "Residuals",
  main = "Residuals vs Fitted"
)
abline(h = 0, lty = 2)

# 2. Residuals versus Thermal Stress
plot(
  dat$ThermalStress,
  residuals_pl,
  pch = 19,
  xlab = "Thermal Stress",
  ylab = "Residuals",
  main = "Residuals vs Thermal Stress"
)
abline(h = 0, lty = 2)

# 3. Distribution of residuals
hist(
  residuals_pl,
  breaks = 12,
  main = "Distribution of Residuals",
  xlab = "Residuals"
)

# 4. Normal Q-Q plot
qqnorm(
  residuals_pl,
  main = "Normal Q-Q Plot"
)
qqline(residuals_pl)

dev.off()
