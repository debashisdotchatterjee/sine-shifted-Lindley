# ---- Libraries ----
library(pracma)
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

R_SGSL <- function(x, theta, mu) {
  1 - G_SGSL(x, theta, mu)
}

h_SGSL <- function(x, theta, mu) {
  g_SGSL(x, theta, mu) / R_SGSL(x, theta, mu)
}


# ---- Parameter Exploration and Calculations ----

theta_values <- c(2, 5)
mu_values <- c(0.5, 1, 2)
x_values <- seq(0, 10, by = 0.01)
colors <- c("blue", "green", "red")

# Create data frame for plotting
results_df <- expand_grid(theta = theta_values, mu = mu_values, x = x_values) %>%
  mutate(
    R = R_SGSL(x, theta, mu),
    h = h_SGSL(x, theta, mu)
  ) %>%
  group_by(theta, mu) %>%
  mutate(
    MTTF = integrate(function(x) x * g_SGSL(x, theta[1], mu[1]), mu[1], Inf)$value
  )

# ---- Plotting ----

# Reliability Function Plot
ggplot(results_df, aes(x = x, y = R, color = factor(mu))) +
  geom_line(size = 1) +
  facet_wrap(~ theta, ncol = 2, scales = "free", labeller = label_bquote(theta==.(theta))) +
  labs(x = "Time (x)", y = "Reliability Function", color = "mu") +
  ggtitle("Reliability Functions of SGSL Distribution") +
  theme_minimal()

# Hazard Function Plot
ggplot(results_df, aes(x = x, y = h, color = factor(mu))) +
  geom_line(size = 1) +
  facet_wrap(~ theta, ncol = 2, scales = "free", labeller = label_bquote(theta==.(theta))) +
  labs(x = "Time (x)", y = "Hazard Function", color = "mu") +
  ggtitle("Hazard Functions of SGSL Distribution") +
  theme_minimal()

# MTTF Bar Plot
results_df %>%
  select(theta, mu, MTTF) %>%
  distinct() %>%
  ggplot(aes(x = factor(mu), y = MTTF, fill = factor(theta))) +
  geom_col(position = "dodge") +
  labs(x = "mu", y = "Mean Time to Failure (MTTF)", fill = "theta") +
  ggtitle("Mean Time to Failure (MTTF) of SGSL Distribution") +
  theme_minimal()
