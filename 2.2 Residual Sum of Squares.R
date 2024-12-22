# Load necessary libraries
library(MASS)

# Load your datasets
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine time and gene expression data
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Assuming X2 is Y output

# Define a function to compute Residual Sum of Squares (RSS)
computeRSS <- function(y, X, beta) {
  y_hat <- X %*% beta  # Predicted values
  rss <- sum((y - y_hat)^2)  # Sum of squared residuals
  return(rss)
}

# Define function to estimate theta using Least Squares
estimate_theta <- function(X, y) {
  theta <- solve(t(X) %*% X) %*% t(X) %*% y
  return(theta)
}

# --- Model 1: y ~ X4 + X3^2 ---
X1 <- as.matrix(cbind(1, data$X4, data$X3^2))  # Design matrix with intercept
y <- as.matrix(data$Y)  # Output variable
theta1 <- estimate_theta(X1, y)  # Estimate parameters
rss1 <- computeRSS(y, X1, theta1)  # Compute RSS
cat("Model 1 RSS:", rss1, "\n")

# --- Model 2: y ~ X4 + X3^2 + X5 ---
X2 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5))
theta2 <- estimate_theta(X2, y)
rss2 <- computeRSS(y, X2, theta2)
cat("Model 2 RSS:", rss2, "\n")

# --- Model 3: y ~ X3 + X4 + X5^3 ---
X3 <- as.matrix(cbind(1, data$X3, data$X4, data$X5^3))
theta3 <- estimate_theta(X3, y)
rss3 <- computeRSS(y, X3, theta3)
cat("Model 3 RSS:", rss3, "\n")

# --- Model 4: y ~ X4 + X3^2 + X5^3 ---
X4 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5^3))
theta4 <- estimate_theta(X4, y)
rss4 <- computeRSS(y, X4, theta4)
cat("Model 4 RSS:", rss4, "\n")

# --- Model 5: y ~ X4 + X1^2 + X3^2 ---
X5 <- as.matrix(cbind(1, data$X4, data$X1^2, data$X3^2))
theta5 <- estimate_theta(X5, y)
rss5 <- computeRSS(y, X5, theta5)
cat("Model 5 RSS:", rss5, "\n")

# Print summary
#cat("\nSummary of Residual Sum of Squares (RSS):\n")
#cat("Model 1 RSS:", rss1, "\n")
#cat("Model 2 RSS:", rss2, "\n")
#cat("Model 3 RSS:", rss3, "\n")
#cat("Model 4 RSS:", rss4, "\n")
#cat("Model 5 RSS:", rss5, "\n")
