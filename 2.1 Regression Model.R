# Load the dataset
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine time and gene expression data
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Rename columns

# Define function to calculate Least Squares parameters
estimate_theta <- function(X, y) {
  theta <- solve(t(X) %*% X) %*% t(X) %*% y  # Compute (X^T X)^-1 X^T y
  return(theta)
}

# Model 1: y ~ X4 + X3^2
X1 <- as.matrix(cbind(1, data$X4, data$X3^2))  # Add bias term (Intercept)
y <- as.matrix(data$Y)
theta1 <- estimate_theta(X1, y)
cat("Model 1 Parameters:\n")
print(theta1)

# Model 2: y ~ X4 + X3^2 + X5
X2 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5))
theta2 <- estimate_theta(X2, y)
cat("\nModel 2 Parameters:\n")
print(theta2)

# Model 3: y ~ X3 + X4 + X5^3
X3 <- as.matrix(cbind(1, data$X3, data$X4, data$X5^3))
theta3 <- estimate_theta(X3, y)
cat("\nModel 3 Parameters:\n")
print(theta3)

# Model 4: y ~ X4 + X3^2 + X5^3
X4 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5^3))
theta4 <- estimate_theta(X4, y)
cat("\nModel 4 Parameters:\n")
print(theta4)

# Model 5: y ~ X4 + X1^2 + X3^2
X5 <- as.matrix(cbind(1, data$X4, data$X1^2, data$X3^2))
theta5 <- estimate_theta(X5, y)
cat("\nModel 5 Parameters:\n")
print(theta5)
