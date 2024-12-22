# Load necessary libraries
library(ggplot2)
library(MASS)

# Load the dataset
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine datasets
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Assuming X2 is Y

# Define Training Output (Y) and Input (X) for Selected Model 4: y ~ X4 + X3^2 + X5^3
y <- as.matrix(data$Y)
X <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5^3))

# Step 1: Least Squares Estimation of Model Parameters ----------------------
theta <- solve(t(X) %*% X) %*% t(X) %*% y  # Estimate parameters
colnames(theta) <- "Value"
cat("Estimated Parameters (Model 4):\n")
print(theta)

# Identify the two parameters with the largest absolute values
param_indices <- order(abs(theta[-1]), decreasing = TRUE)[1:2] + 1  # Exclude Intercept
param_names <- colnames(X)[param_indices]
cat("Parameters selected for ABC:", param_names, "\n")

# Fix the other parameters as constants
theta_fixed <- theta
theta_fixed[-param_indices] <- theta[-param_indices]  # Fix parameters except selected ones

# Step 2: Define Uniform Priors for Selected Parameters ---------------------
set.seed(123)
n_samples <- 10000  # Number of prior samples

# Define prior ranges (± 20% around estimated values)
prior_ranges <- list(
  param1 = c(theta[param_indices[1]] * 0.8, theta[param_indices[1]] * 1.2),
  param2 = c(theta[param_indices[2]] * 0.8, theta[param_indices[2]] * 1.2)
)

# Draw samples from Uniform priors
param1_samples <- runif(n_samples, min = prior_ranges$param1[1], max = prior_ranges$param1[2])
param2_samples <- runif(n_samples, min = prior_ranges$param2[1], max = prior_ranges$param2[2])

# Step 3: Perform Rejection ABC ---------------------------------------------
# Tolerance level for ABC
tolerance <- 0.01

# Function to compute RSS for given parameter values
computeRSS <- function(y, X, param1, param2, fixed_params) {
  beta <- fixed_params
  beta[param_indices[1]] <- param1
  beta[param_indices[2]] <- param2
  y_hat <- X %*% beta
  rss <- sum((y - y_hat)^2)
  return(rss)
}

# Initialize acceptance storage
accepted_param1 <- c()
accepted_param2 <- c()

# Perform rejection ABC
for (i in 1:n_samples) {
  rss <- computeRSS(y, X, param1_samples[i], param2_samples[i], theta_fixed)
  if (rss < tolerance) {
    accepted_param1 <- c(accepted_param1, param1_samples[i])
    accepted_param2 <- c(accepted_param2, param2_samples[i])
  }
}

cat("Number of accepted samples:", length(accepted_param1), "\n")

# Step 4: Plot Joint and Marginal Posterior Distributions -------------------
# Combine accepted samples
accepted_data <- data.frame(param1 = accepted_param1, param2 = accepted_param2)

# Joint Posterior Distribution
joint_plot <- ggplot(accepted_data, aes(x = param1, y = param2)) +
  geom_point(alpha = 0.5, color = "blue") +
  labs(title = "Joint Posterior Distribution", x = param_names[1], y = param_names[2]) +
  theme_minimal()

# Marginal Posterior for Parameter 1
marginal1_plot <- ggplot(accepted_data, aes(x = param1)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = paste("Marginal Posterior of", param_names[1]), x = param_names[1], y = "Density") +
  theme_minimal()

# Marginal Posterior for Parameter 2
marginal2_plot <- ggplot(accepted_data, aes(x = param2)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  labs(title = paste("Marginal Posterior of", param_names[2]), x = param_names[2], y = "Density") +
  theme_minimal()

# Display Plots
gridExtra::grid.arrange(joint_plot, marginal1_plot, marginal2_plot, nrow = 2)

# Save Plots
ggsave("E:\\Rpractice\\regression\\joint_posterior.png", joint_plot, width = 6, height = 6)
ggsave("E:\\Rpractice\\regression\\marginal1_posterior.png", marginal1_plot, width = 6, height = 4)
ggsave("E:\\Rpractice\\regression\\marginal2_posterior.png", marginal2_plot, width = 6, height = 4)
