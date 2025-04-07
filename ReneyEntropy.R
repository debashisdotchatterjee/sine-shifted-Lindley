# Load necessary library
library(pracma)

# Define the shifted Lindley CDF
F_shifted_lindley <- function(x, theta, mu) {
  1 - ((1 + theta * (1 + x)) / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

# Define the shifted Lindley PDF
f_shifted_lindley <- function(x, theta, mu) {
  (theta^2 / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

# Define the SGSL PDF
g_SGSL <- function(x, theta, mu) {
  (pi / 2) * cos((pi / 2) * F_shifted_lindley(x, theta, mu)) * f_shifted_lindley(x, theta, mu)
}

# Define the integrand for the R\'enyi entropy
renyi_integrand <- function(x, alpha, theta, mu) {
  ((pi / 2) * cos((pi / 2) * F_shifted_lindley(x, theta, mu)) * f_shifted_lindley(x, theta, mu))^alpha
}

# Define the R\'enyi entropy function
H_alpha_SGSL <- function(alpha, theta, mu) {
  integral_value <- integral(renyi_integrand, mu, Inf, alpha = alpha, theta = theta, mu = mu)
  H_alpha <- (1 / (1 - alpha)) * log((pi / 2)^alpha * (theta^2 / (1 + theta * (1 + mu)))^alpha * integral_value)
  return(H_alpha)
}

# Example parameters
alpha <- 2
theta <- 2
mu <- 1

# Calculate the R\'enyi entropy
H_alpha_value <- H_alpha_SGSL(alpha, theta, mu)
print(H_alpha_value)

#####################################################

# ---- Libraries ----
library(pracma)
library(ggplot2)
library(tidyverse)

# ---- Define Functions (with support restriction) ----
F_shifted_lindley <- function(x, theta, mu) {
  ifelse(x > mu, 1 - ((1 + theta * (1 + x)) / (1 + theta * (1 + mu))) * exp(-theta * (x - mu)), 0)
}

f_shifted_lindley <- function(x, theta, mu) {
  ifelse(x > mu, (theta^2 / (1 + theta * (1 + mu))) * exp(-theta * (x - mu)), 0)
}

G_SGSL <- function(x, theta, mu) {
  sin((pi / 2) * F_shifted_lindley(x, theta, mu))
}

g_SGSL <- function(x, theta, mu) {
  ifelse(x > mu, (pi / 2) * abs(cos((pi / 2) * F_shifted_lindley(x, theta, mu))) * f_shifted_lindley(x, theta, mu), 0)
}

# ---- Rényi Entropy Integrand ----
renyi_integrand <- function(x, alpha, theta, mu) {
  g_SGSL(x, theta, mu)^alpha
}

# ---- Rényi Entropy Function (Vectorized) ----
H_alpha_SGSL <- function(alpha_values, theta, mu) {
  sapply(alpha_values, function(alpha) {
    integral_value <- integrate(renyi_integrand, mu, Inf, alpha = alpha, theta = theta, mu = mu)$value
    (1 / (1 - alpha)) * log((pi / 2)^alpha * (theta^2 / (1 + theta * (1 + mu)))^alpha * integral_value)
  })
}

# ---- Parameter Values ----
theta_values <- c(2, 5)  # Different theta values for comparison
mu_values <- c(0.5, 1, 2) # Different mu values for comparison
alpha_values <- seq(0.1, 10, by = 0.1)  # Range of alpha values

# ---- Calculate Rényi Entropy ----
renyi_results <- expand_grid(theta = theta_values, mu = mu_values, alpha = alpha_values) %>%
  rowwise() %>%
  mutate(renyi_entropy = H_alpha_SGSL(alpha, theta, mu)) %>%
  ungroup()

# ---- Plotting ----

# Base plot with facets for theta and mu
base_plot <- ggplot(renyi_results, aes(x = alpha, y = renyi_entropy, color = factor(mu))) +
  geom_line(size = 1) +
  labs(x = expression(alpha), y = "Rényi Entropy", color = expression(mu)) +
  theme_minimal() +
  facet_grid(theta ~ mu, labeller = label_bquote(theta==.(theta), mu==.(mu)))

# Add vertical line at alpha = 1 for reference (Shannon Entropy)
final_plot <- base_plot + geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  ggtitle("Rényi Entropy of SGSL Distribution for Varying Parameters")

print(final_plot)
