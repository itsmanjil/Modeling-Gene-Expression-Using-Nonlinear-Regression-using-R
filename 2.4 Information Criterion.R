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

computeLogLikelihood <- function(n, rss) {
  sigma2 <- rss / (n - 1)  # Variance estimate
  log_likelihood <- -(n / 2) * log(2 * pi) - (n / 2) * log(sigma2) - (rss / (2 * sigma2))
  return(log_likelihood)
}

# Function to compute AIC and BIC
computeAIC_BIC <- function(log_likelihood, k, n) {
  AIC <- 2 * k - 2 * log_likelihood
  BIC <- k * log(n) - 2 * log_likelihood
  return(list(AIC = AIC, BIC = BIC))
}

# Initialize variables
y <- as.matrix(data$Y)  # Output variable
n <- nrow(data)         # Number of samples

# --- Model 1: y ~ X4 + X3^2 ---
X1 <- as.matrix(cbind(1, data$X4, data$X3^2))  # Design matrix
theta1 <- estimate_theta(X1, y)
rss1 <- computeRSS(y, X1, theta1)
logLik1 <- computeLogLikelihood(n, rss1)
metrics1 <- computeAIC_BIC(logLik1, ncol(X1), n)
cat("Model 1 AIC:", metrics1$AIC, "BIC:", metrics1$BIC, "\n")

# --- Model 2: y ~ X4 + X3^2 + X5 ---
X2 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5))
theta2 <- estimate_theta(X2, y)
rss2 <- computeRSS(y, X2, theta2)
logLik2 <- computeLogLikelihood(n, rss2)
metrics2 <- computeAIC_BIC(logLik2, ncol(X2), n)
cat("Model 2 AIC:", metrics2$AIC, "BIC:", metrics2$BIC, "\n")

# --- Model 3: y ~ X3 + X4 + X5^3 ---
X3 <- as.matrix(cbind(1, data$X3, data$X4, data$X5^3))
theta3 <- estimate_theta(X3, y)
rss3 <- computeRSS(y, X3, theta3)
logLik3 <- computeLogLikelihood(n, rss3)
metrics3 <- computeAIC_BIC(logLik3, ncol(X3), n)
cat("Model 3 AIC:", metrics3$AIC, "BIC:", metrics3$BIC, "\n")

# --- Model 4: y ~ X4 + X3^2 + X5^3 ---
X4 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5^3))
theta4 <- estimate_theta(X4, y)
rss4 <- computeRSS(y, X4, theta4)
logLik4 <- computeLogLikelihood(n, rss4)
metrics4 <- computeAIC_BIC(logLik4, ncol(X4), n)
cat("Model 4 AIC:", metrics4$AIC, "BIC:", metrics4$BIC, "\n")

# --- Model 5: y ~ X4 + X1^2 + X3^2 ---
X5 <- as.matrix(cbind(1, data$X4, data$X1^2, data$X3^2))
theta5 <- estimate_theta(X5, y)
rss5 <- computeRSS(y, X5, theta5)
logLik5 <- computeLogLikelihood(n, rss5)
metrics5 <- computeAIC_BIC(logLik5, ncol(X5), n)
cat("Model 5 AIC:", metrics5$AIC, "BIC:", metrics5$BIC, "\n")

# Summary of AIC and BIC
cat("\nSummary of AIC and BIC:\n")
cat("Model 1 AIC:", metrics1$AIC, "BIC:", metrics1$BIC, "\n")
cat("Model 2 AIC:", metrics2$AIC, "BIC:", metrics2$BIC, "\n")
cat("Model 3 AIC:", metrics3$AIC, "BIC:", metrics3$BIC, "\n")
cat("Model 4 AIC:", metrics4$AIC, "BIC:", metrics4$BIC, "\n")
cat("Model 5 AIC:", metrics5$AIC, "BIC:", metrics5$BIC, "\n")


# Combine results into a data frame
model_results <- data.frame(
  Models = c("Model 1", "Model 2", "Model 3", "Model 4", "Model 5"),
  BIC = c(metrics1$BIC, metrics2$BIC, metrics3$BIC, metrics4$BIC, metrics5$BIC)
)

# Print the data frame
print(model_results)

