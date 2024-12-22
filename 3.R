# Load necessary libraries
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)
if (!require(gridExtra)) install.packages("gridExtra", dependencies = TRUE)

# Load the dataset
gene_data <- read.csv("E:/Rpractice/regression/data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:/Rpractice/regression/time_1673241270748.csv")

# Combine time and gene expression data
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")

# Assume `result` contains parameter estimates from Task 2
result <- c(-15.5, 9.8, 1.2, -3.4, 0.8)  # Replace these with the actual values

# Design matrix for the regression model
X <- as.matrix(cbind(1, data$X4, data$X1^2, data$X3^2))  # Include intercept and predictors
observed_data <- as.matrix(data$Y)  # Response variable (observed Y values)

# Debugging information
cat("Design matrix dimensions (rows x columns):", dim(X), "\n")

# Select the two parameters with the largest absolute values
param_indices <- order(abs(result), decreasing = TRUE)[1:2]
estimated_params <- result[param_indices]  # Selected parameters for ABC
fixed_params <- result[-param_indices]    # Remaining fixed parameters

# Validate the length of fixed_params dynamically
num_predictors <- ncol(X)  # Number of predictors in the design matrix
if (length(fixed_params) != (num_predictors - 2)) {
  # Adjust fixed_params to match the expected length dynamically
  fixed_params <- result[-param_indices][1:(num_predictors - 2)]
}

# Print parameters for debugging
cat("Selected parameters for ABC (estimated_params):", estimated_params, "\n")
cat("Fixed parameters (fixed_params):", fixed_params, "\n")

# Define prior distributions for the selected parameters
set.seed(123)
num_samples <- 1000
prior_theta1 <- runif(num_samples, min = estimated_params[1] - 5, max = estimated_params[1] + 5)
prior_theta2 <- runif(num_samples, min = estimated_params[2] - 5, max = estimated_params[2] + 5)

# Rejection ABC process
epsilon <- 2  # Tolerance threshold for RSS
accepted_theta1 <- c()
accepted_theta2 <- c()

for (i in 1:num_samples) {
  theta <- c(prior_theta1[i], prior_theta2[i], fixed_params)
  simulated_data <- X %*% matrix(theta, ncol = 1) + rnorm(nrow(X), mean = 0, sd = 1)
  rss <- sum((observed_data - simulated_data)^2)
  
  if (rss < epsilon) {
    accepted_theta1 <- c(accepted_theta1, prior_theta1[i])
    accepted_theta2 <- c(accepted_theta2, prior_theta2[i])
  }
}

# Check if any samples were accepted
if (length(accepted_theta1) == 0 || length(accepted_theta2) == 0) {
  stop("No samples accepted. Adjust the epsilon value or priors.")
}

# Create a data frame for accepted values
accepted_data <- data.frame(theta1 = accepted_theta1, theta2 = accepted_theta2)

# Joint Posterior Distribution Plot
joint_plot <- ggplot(accepted_data, aes(x = theta1, y = theta2)) +
  geom_point(alpha = 0.6, color = "blue") +
  labs(title = "Joint and Marginal Posterior Distribution", x = "Theta 1 (thetabias)", y = "Theta 2 (theta X4)") +
  theme_minimal()

# Marginal Posterior Distributions
marginal_theta1_plot <- ggplot(accepted_data, aes(x = theta1)) +
  geom_density(fill = "blue", alpha = 0.4) +
  labs(title = "Marginal Posterior Distribution of Theta 1", x = "Theta 1 (thetabias)", y = "Density") +
  theme_minimal()

marginal_theta2_plot <- ggplot(accepted_data, aes(x = theta2)) +
  geom_density(fill = "green", alpha = 0.4) +
  labs(title = "Marginal Posterior Distribution of Theta 2", x = "Theta 2 (theta X4)", y = "Density") +
  theme_minimal()

# Save plots as an image
save_path <- "posterior_distributions.png"
g <- grid.arrange(joint_plot, marginal_theta1_plot, marginal_theta2_plot, nrow = 2)
ggsave(filename = save_path, plot = g, width = 10, height = 8, dpi = 300)

# Print success message
cat("Posterior distributions saved to:", save_path, "\n")

# Display plots
print(joint_plot)
print(marginal_theta1_plot)
print(marginal_theta2_plot)
