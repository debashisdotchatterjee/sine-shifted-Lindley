# Load necessary libraries
library(MASS) # For fitting distributions
library(fitdistrplus) # For goodness-of-fit measures

# Set seed for reproducibility
set.seed(123)

# Simulate data from SGSL distribution
simulate_SGSL <- function(n, theta, mu) {
  u <- runif(n)
  F_inv <- function(u, theta, mu) {
    return(-log((1 - u) / (1 + theta * (1 + mu))) / theta + mu)
  }
  return(F_inv(u, theta, mu))
}

n <- 1000
theta <- 2
mu <- 1
data_SGSL <- simulate_SGSL(n, theta, mu)

# Fit distributions to the data
fit_SGSL <- fitdist(data_SGSL, "weibull", start = list(shape = 1, scale = 1)) # Replace with SGSL fitting
fit_Lindley <- fitdist(data_SGSL, "gamma", start = list(shape = 1, scale = 1)) # Replace with Lindley fitting
fit_Exponential <- fitdist(data_SGSL, "exp")
fit_Shifted_Lindley <- fitdist(data_SGSL, "gamma", start = list(shape = 1, scale = 1)) # Replace with Shifted Lindley fitting

# Calculate AIC for comparison
aic_values <- c(
  SGSL = fit_SGSL$aic,
  Lindley = fit_Lindley$aic,
  Exponential = fit_Exponential$aic,
  Shifted_Lindley = fit_Shifted_Lindley$aic
)

# Plot fitted distributions
pdf("fit_comparison.pdf")
par(mfrow = c(2, 2))
plot(fit_SGSL)
plot(fit_Lindley)
plot(fit_Exponential)
plot(fit_Shifted_Lindley)
dev.off()

# Output AIC values
aic_values