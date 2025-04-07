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

# Define the MGF integrand
mgf_integrand <- function(x, t, theta, mu) {
  e^(t * x) * g_SGSL(x, theta, mu)
}

# Define the MGF function
M_SGSL <- function(t, theta, mu) {
  integral(mgf_integrand, mu, Inf, t = t, theta = theta, mu = mu)
}

# Define the quantile function (solve for x)
quantile_function <- function(p, theta, mu) {
  y <- 1 - 2 / pi * asin(p)
  equation <- function(x) {
    log(1 + theta * (1 + x)) - log(1 + theta * (1 + mu)) - log(y) - theta * (x - mu)
  }
  uniroot(equation, lower = mu, upper = mu + 100)$root
}

# Example parameters
theta <- 2
mu <- 1
t <- 0.5
p <- 0.5

# Calculate the MGF value
mgf_value <- M_SGSL(t, theta, mu)
print(mgf_value)

# Calculate the quantile function value
quantile_value <- quantile_function(p, theta, mu)
print(quantile_value)
################################################

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

# ---- Reliability and Hazard Functions ----
R_SGSL <- function(x, theta, mu) {
  1 - G_SGSL(x, theta, mu)
}

h_SGSL <- function(x, theta, mu) {
  g_SGSL(x, theta, mu) / R_SGSL(x, theta, mu)
}

# ---- MGF Integrand ----
mgf_integrand <- function(x, t, theta, mu) {
  exp(t * x) * g_SGSL(x, theta, mu)
}

# ---- MGF Function (Monte Carlo Approximation) ----
M_SGSL_MC <- function(t, theta, mu, n_simulations = 10000) {
  samples <- rSGSL(n_simulations, theta, mu)
  mean(exp(t * samples))
}

# ---- Quantile Function (Solve for x) ----
Q_SGSL <- function(p, theta, mu) {
  y <- 1 - (2 / pi) * asin(p)
  equation <- function(x) F_shifted_lindley(x, theta, mu) - (2 / pi) * asin(p)
  uniroot(equation, c(mu, mu + 100), tol = .Machine$double.eps^0.5)$root
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
theta_values <- seq(0.5, 5, by = 0.5)
mu_values <- seq(0, 2, by = 0.5)
t_values <- seq(-1, 1, by = 0.1)
p_values <- seq(0, 1, by = 0.01)
alpha_values <- seq(0.1, 10, by = 0.1)

# ---- Calculate Values for Plots ----

# MGF
mgf_results <- expand_grid(theta = theta_values, mu = mu_values, t = t_values) %>%
  rowwise() %>%
  mutate(mgf_value = M_SGSL_MC(t, theta, mu)) 

# Quantiles
quantile_results <- expand_grid(theta = theta_values, mu = mu_values, p = p_values) %>%
  rowwise() %>%
  mutate(quantile = Q_SGSL(p, theta, mu))

# Reliability and Hazard
x_values <- seq(0, 10, by = 0.01) 
rel_haz_results <- expand_grid(theta = theta_values, mu = mu_values, x = x_values) %>%
  mutate(R = R_SGSL(x, theta, mu), h = h_SGSL(x, theta, mu))

# Renyi Entropy
renyi_results <- expand_grid(theta = theta_values, mu = mu_values, alpha = alpha_values) %>%
  rowwise() %>%
  mutate(renyi_entropy = H_alpha_SGSL(alpha, theta, mu)) %>%
  ungroup()


# ---- Plotting ----

# MGF Plot
ggplot(mgf_results, aes(x = t, y = mgf_value, color = factor(mu))) +
  geom_line(size = 1) +
  facet_wrap(~ theta, ncol = 2, labeller = label_bquote(theta==.(theta))) +
  labs(x = "t", y = "MGF", color = expression(mu)) +
  ggtitle("Moment Generating Function (MGF) of SGSL") +
  theme_minimal()

# Quantile Plot
ggplot(quantile_results, aes(x = p, y = quantile, color = factor(mu))) +
  geom_line(size = 1) +
  facet_wrap(~ theta, ncol = 2, labeller = label_bquote(theta==.(theta))) +
  labs(x = "p", y = "Quantile", color = expression(mu)) +
  ggtitle("Quantile Function of SGSL") +
  theme_minimal()

# Reliability Function Plot
ggplot(rel_haz_results, aes(x = x, y = R, color = factor(mu))) +
  geom_line(size = 1) +
  facet_wrap(~ theta, ncol = 2, scales = "free", labeller = label_bquote(theta==.(theta))) +
  labs(x = "Time (x)", y = "Reliability Function", color = "mu") +
  ggtitle("Reliability Functions of SGSL Distribution") +
  theme_minimal()

# Hazard Function Plot
ggplot(rel_haz_results, aes(x = x, y = h, color = factor(mu))) +
  geom_line(size = 1) +
  facet_wrap(~ theta, ncol = 2, scales = "free", labeller = label_bquote(theta==.(theta))) +
  labs(x = "Time (x)", y = "Hazard Function", color = "mu") +
  ggtitle("Hazard Functions of SGSL Distribution") +
  theme_minimal()

# Rényi Entropy Plot
base_plot <- ggplot(renyi_results, aes(x = alpha, y = renyi_entropy, color = factor(mu))) +
  geom_line(size = 1) +
  labs(x = expression(alpha), y = "Rényi Entropy", color = expression(mu)) +
  theme_minimal() +
  facet_grid(theta ~ mu, labeller = label_bquote(theta==.(theta), mu==.(mu)))

# Add vertical line at alpha = 1 for reference (Shannon Entropy)
final_plot <- base_plot + geom_vline(xintercept = 1, linetype = "dashed", color = "gray") +
  ggtitle("Rényi Entropy of SGSL Distribution for Varying Parameters")

print(final_plot)

# Heatmaps for MGF and Quantile Functions
mgf_heatmap_data <- mgf_results %>%
  group_by(theta, mu) %>%
  summarize(mean_mgf = mean(mgf_value))

quantile_heatmap_data <- quantile_results %>%
  group_by(theta, mu) %>%
  summarize(mean_quantile = mean(quantile))

# Plot Heatmaps
ggplot(mgf_heatmap_data, aes(x = mu, y = theta, fill = mean_mgf)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
  labs(x = expression(mu), y = expression(theta), fill = "Mean MGF") +
  ggtitle("Heatmap of Mean MGF for SGSL")

ggplot(quantile_heatmap_data, aes(x = mu,
                                  

