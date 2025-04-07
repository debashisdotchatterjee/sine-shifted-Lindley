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

# ---- Inverse CDF of Shifted Lindley (Using uniroot for numerical inversion) ----

inverse_F_shifted_lindley <- function(p, theta, mu) {
  f <- function(x) F_shifted_lindley(x, theta, mu) - p
  uniroot(f, c(mu, mu + 100), tol = .Machine$double.eps^0.5)$root  # Start search from mu
}

# ---- Inverse Transform Sampling Function ----

rSGSL <- function(n, theta, mu) {
  U <- runif(n)  # Generate uniform random variables
  sapply(U, function(u) inverse_F_shifted_lindley((2/pi) * asin(u), theta, mu))
}


# ---- Simulation & Plotting ----

theta <- 3
mu_values <- c(0, 0.5, 1) 
n <- 10000  # Number of random samples

# Plotting
par(mfrow = c(1, 2))  # Arrange plots in a row

for (mu in mu_values) {
  # Generate samples
  samples <- rSGSL(n, theta, mu)

  # Calculate theoretical PDF values for plotting
  x_values <- seq(mu, 10, by = 0.01)  # Start x_values from mu
  pdf_values <- g_SGSL(x_values, theta, mu)

  # Histogram of Simulated Samples
  hist(samples, breaks = 30, freq = FALSE, col = "skyblue",
       xlab = "x", ylab = "Density", main = paste("SGSL Samples (theta = 3, mu =", mu, ")"))
  lines(x_values, pdf_values, col = "red", lwd = 2)  # Overlay theoretical PDF

  # Empirical CDF vs. Theoretical CDF
  plot(ecdf(samples), col = "blue", pch = ".", verticals = TRUE, do.points = FALSE,
       xlab = "x", ylab = "Cumulative Probability", main = paste("SGSL CDF (theta = 3, mu =", mu, ")"))
  lines(x_values, G_SGSL(x_values, theta, mu), col = "red", lwd = 2)
}

