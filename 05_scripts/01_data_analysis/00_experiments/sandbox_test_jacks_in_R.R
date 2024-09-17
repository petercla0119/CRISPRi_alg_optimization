#######################################
# IGNORE LINES BELOW - USED AS SANDBOX
#######################################

# Test JACKS in R

# Load necessary library
library(dplyr)
library(zoo)

# Function to compute the observed log2 fold change
compute_log2_fold_change <- function(treatment_counts, control_counts) {
  # Add pseudocount of 32 and log-transform the counts
  log_transform <- function(x) log2(x + 32)
  # Normalize by the median across all guides in each replicate
  median_normalize <- function(counts) {
    log_counts <- apply(counts, 2, log_transform)
    medians <- apply(log_counts, 2, median)
    normalized_counts <- sweep(log_counts, 2, medians, '-')
    return(normalized_counts)
  }
  # Normalize treatment and control counts
  T_normalized <- median_normalize(treatment_counts)
  C_normalized <- median_normalize(control_counts)
  # Compute the mean across replicates for treatment and control
  T_mean <- rowMeans(T_normalized)
  C_mean <- rowMeans(C_normalized)
  # Compute the observed log2 fold change
  log2_fold_change <- T_mean - C_mean
  return(log2_fold_change)
}
# Log-transform and median-normalize function ####
log_transform_and_normalize <- function(counts) {
  log_counts <- log2(counts + 32)
  medians <- apply(log_counts, 2, median)
  normalized_counts <- sweep(log_counts, 2, medians, '-')
  return(normalized_counts)
}
## Apply normalization ####
T_normalized <- log_transform_and_normalize(treatment_counts)
C_normalized <- log_transform_and_normalize(control_counts)

# Compute mean and variance function ####
compute_mean_variance <- function(normalized_counts) {
  means <- rowMeans(normalized_counts)
  variances <- apply(normalized_counts, 1, var)
  return(list(means = means, variances = variances))
}
# Compute mean and variance ####
T_stats <- compute_mean_variance(T_normalized)
C_stats <- compute_mean_variance(C_normalized)

# Combine treatment and control statistics ####
mean_variance_pairs <- data.frame(
  mean_value = c(T_stats$means, C_stats$means),
  variance = c(T_stats$variances, C_stats$variances)
)
# Sort by mean value
sorted_mean_variance <- mean_variance_pairs %>% arrange(mean_value)
# Apply moving average filter
apply_moving_average <- function(sorted_variances, window_size = 800) {
  filter_size <- min(window_size, length(sorted_variances))
  smoothed_variances <- zoo::rollmean(sorted_variances, k = filter_size, fill = NA, align = "center")
  smoothed_variances <- na.approx(smoothed_variances, method = "constant", rule = 2)
  return(smoothed_variances)
}
# Assuming window size of 800
sorted_mean_variance$smoothed_variance <- apply_moving_average(sorted_mean_variance$variance)


# Placeholder values for the priors (replace with actual values as needed)
mu_x <- 0
sigma2_x <- 1
mu_w <- 0
sigma2_w <- 1000
kappa <- 0.5
# Calculate the parameters for the Gamma distribution (τ_g,i,l)
compute_tau_parameters <- function(variances, kappa) {
  tau_parameters <- kappa * variances
  return(tau_parameters)
}
# Compute τ_g,i,l parameters for treatment and control
tau_T <- compute_tau_parameters(T_stats$variances, kappa)
tau_C <- compute_tau_parameters(C_stats$variances, kappa)
# Define the function to compute β and τ* parameters
compute_beta_tau_star <- function(y_gil, x_gi, w_gil, tau_gil) {
  beta_star <- y_gil^2 - 2 * y_gil * x_gi * w_gil + x_gi^2 * w_gil^2
  tau_star <- 0.5 * sum(y_gil^2) + tau_gil * beta_star
  return(list(beta_star = beta_star, tau_star = tau_star))
}
# Placeholder values for y_gil, x_gi, w_gil, tau_gil (replace with actual values)
y_gil <- c(1, 2, 3)  # Example values
x_gi <- 1
w_gil <- 1
tau_gil <- 1
# Compute β* and τ*
beta_tau_star <- compute_beta_tau_star(y_gil, x_gi, w_gil, tau_gil)
# Compute posterior distributions for x_gi, w_gi,l, and τ_g,i,l (Q_x, Q_w, Q_tau)
compute_posterior_distributions <- function(mu, sigma2, y_gil, x_gi, w_gil, tau_gil, beta_star) {
  Q_x_gi <- dnorm(x_gi, mean = mu, sd = sqrt(sigma2))
  Q_w_gil <- dnorm(w_gil, mean = mu, sd = sqrt(sigma2))
  Q_tau_gil <- dgamma(tau_gil, shape = 0.5, rate = 0.5 * beta_star)
  return(list(Q_x_gi = Q_x_gi, Q_w_gil = Q_w_gil, Q_tau_gil = Q_tau_gil))
}
# Compute posterior distributions
posterior_distributions <- compute_posterior_distributions(mu_x, sigma2_x, y_gil, x_gi, w_gil, tau_gil, beta_tau_star$beta_star)
# Print results (example)
print(posterior_distributions)


# Test Data ####

# Example data
treatment_counts <- matrix(c(
  1000, 1200, 1100,
  1500, 1300, 1400,
  900,  1100, 950,
  1200, 1250, 1150
), nrow = 4, byrow = TRUE)
colnames(treatment_counts) <- c("replicate1", "replicate2", "replicate3")
rownames(treatment_counts) <- c("sgRNA1", "sgRNA2", "sgRNA3", "sgRNA4")

control_counts <- matrix(c(
  1100, 1300, 1200,
  1400, 1350, 1250,
  950,  1050, 1000,
  1150, 1220, 1180
), nrow = 4, byrow = TRUE)
colnames(control_counts) <- c("replicate1", "replicate2", "replicate3")
rownames(control_counts) <- c("sgRNA1", "sgRNA2", "sgRNA3", "sgRNA4")

# Compute the observed log2 fold change
log2_fold_change <- compute_log2_fold_change(treatment_counts, control_counts)
print(log2_fold_change)
