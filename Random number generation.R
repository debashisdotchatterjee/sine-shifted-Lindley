# ---- Define Functions (Same as before) ----

# ---- Define Functions ----

F_shifted_lindley <- function(x, theta, mu) {
  1 - ((1 + theta * (1 + x)) / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

f_shifted_lindley <- function(x, theta, mu) {
  (theta^2 / (1 + theta * (1 + mu))) * exp(-theta * (x - mu))
}

G_SGSL <- function(x, theta, mu) {
  sin((pi / 2) * F_shifted_lindley(x, theta, mu))
}

g_SGSL <- function(x, theta, mu) {
  (pi / 2) * cos((pi / 2) * F_shifted_lindley(x, theta, mu)) * f_shifted_lindley(x, theta, mu)
}

h_SGSL <- function(x, theta, mu) {
  g_SGSL(x, theta, mu) / (1 - G_SGSL(x, theta, mu))
}


#########################

# ---- Inverse CDF of Shifted Lindley (Using uniroot for numerical inversion) ----

inverse_F_shifted_lindley <- function(p, theta, mu) {
  f <- function(x) F_shifted_lindley(x, theta, mu) - p
  uniroot(f, c(mu, mu + 100), tol = .Machine$double.eps^0.5)$root  
}

# ---- Inverse Transform Sampling Function ----

rSGSL <- function(n, theta, mu) {
  U <- runif(n)  # Generate uniform random variables
  sapply(U, function(u) inverse_F_shifted_lindley((2/pi) * asin(u), theta, mu))
}


# ---- Simulation & Plotting ----

theta <- 2
mu <- 1
n <- 10000  # Number of random samples

# Generate samples
samples <- rSGSL(n, theta, mu)

# Calculate theoretical PDF values for plotting
x_values <- seq(mu, 10, by = 0.01)
pdf_values <- g_SGSL(x_values, theta, mu)

# Plotting
par(mfrow = c(1, 2))

# Histogram of Simulated Samples
hist(samples, breaks = 30, freq = FALSE, col = "skyblue",
     xlab = "x", ylab = "Density", main = "SGSL Simulated Samples")
lines(x_values, pdf_values, col = "red", lwd = 2)  # Overlay theoretical PDF

# Empirical CDF vs. Theoretical CDF
plot(ecdf(samples), col = "blue", pch = ".", verticals = TRUE, do.points = FALSE,
     xlab = "x", ylab = "Cumulative Probability", main = "SGSL CDF: Empirical vs. Theoretical")
lines(x_values, G_SGSL(x_values, theta, mu), col = "red", lwd = 2)
grid()
