# Load necessary libraries
library(MASS)

# Load the dataset
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine time and gene expression data
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Assuming X2 is Y output

# Function to compute RSS
computeRSS <- function(y, X, beta) {
  y_hat <- X %*% beta
  rss <- sum((y - y_hat)^2)
  return(rss)
}

# Function to estimate theta
estimate_theta <- function(X, y) {
  theta <- solve(t(X) %*% X) %*% t(X) %*% y
  return(theta)
}

# Function to compute log-likelihood
computeLogLikelihood <- function(n, rss) {
  sigma2 <- rss / (n - 1)  # Variance estimate
  log_likelihood <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (rss / (2 * sigma2))
  return(log_likelihood)
}

# Initialize variables
y <- as.matrix(data$Y)  # Output variable
n <- nrow(data)         # Number of samples

# --- Model 1: y ~ X4 + X3^2 ---
X1 <- as.matrix(cbind(1, data$X4, data$X3^2))  # Design matrix
theta1 <- estimate_theta(X1, y)               # Estimate parameters
rss1 <- computeRSS(y, X1, theta1)             # Compute RSS
logLik1 <- computeLogLikelihood(n, rss1)      # Compute log-likelihood
cat("Model 1 Log-Likelihood:", logLik1, "\n")

# --- Model 2: y ~ X4 + X3^2 + X5 ---
X2 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5))
theta2 <- estimate_theta(X2, y)
rss2 <- computeRSS(y, X2, theta2)
logLik2 <- computeLogLikelihood(n, rss2)
cat("Model 2 Log-Likelihood:", logLik2, "\n")

# --- Model 3: y ~ X3 + X4 + X5^3 ---
X3 <- as.matrix(cbind(1, data$X3, data$X4, data$X5^3))
theta3 <- estimate_theta(X3, y)
rss3 <- computeRSS(y, X3, theta3)
logLik3 <- computeLogLikelihood(n, rss3)
cat("Model 3 Log-Likelihood:", logLik3, "\n")

# --- Model 4: y ~ X4 + X3^2 + X5^3 ---
X4 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5^3))
theta4 <- estimate_theta(X4, y)
rss4 <- computeRSS(y, X4, theta4)
logLik4 <- computeLogLikelihood(n, rss4)
cat("Model 4 Log-Likelihood:", logLik4, "\n")

# --- Model 5: y ~ X4 + X1^2 + X3^2 ---
X5 <- as.matrix(cbind(1, data$X4, data$X1^2, data$X3^2))
theta5 <- estimate_theta(X5, y)
rss5 <- computeRSS(y, X5, theta5)
logLik5 <- computeLogLikelihood(n, rss5)
cat("Model 5 Log-Likelihood:", logLik5, "\n")

# Summary of log-likelihoods
cat("\nSummary of Log-Likelihoods:\n")
cat("Model 1 Log-Likelihood:", logLik1, "\n")
cat("Model 2 Log-Likelihood:", logLik2, "\n")
cat("Model 3 Log-Likelihood:", logLik3, "\n")
cat("Model 4 Log-Likelihood:", logLik4, "\n")
cat("Model 5 Log-Likelihood:", logLik5, "\n")


# Compute variances
variance1 <- rss1 / (n - ncol(X1))
variance2 <- rss2 / (n - ncol(X2))
variance3 <- rss3 / (n - ncol(X3))
variance4 <- rss4 / (n - ncol(X4))
variance5 <- rss5 / (n - ncol(X5))

# Combine results into a data frame
model_results <- data.frame(
  Models = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"),
  Variance = c(variance1, variance2, variance3, variance4, variance5),
  Likelihood = c(logLik1, logLik2, logLik3, logLik4, logLik5)
)

# Print the data frame
print(model_results)

