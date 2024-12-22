# Load necessary libraries
library(MASS)  # For solving (X^T X)^-1 X^T y
library(ggplot2)  # For Q-Q plots

# Load the dataset
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine the datasets
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Assuming X2 is Y output

# Function to compute residuals
computeResiduals <- function(y, X, beta) {
  y_hat <- X %*% beta  # Predicted values
  residuals <- y - y_hat  # Residuals
  return(residuals)
}

# Function to estimate theta
estimate_theta <- function(X, y) {
  theta <- solve(t(X) %*% X) %*% t(X) %*% y
  return(theta)
}

# Initialize variables
y <- as.matrix(data$Y)

# --- Model 1: y ~ X4 + X3^2 ---
X1 <- as.matrix(cbind(1, data$X4, data$X3^2))  # Design matrix
theta1 <- estimate_theta(X1, y)
residuals1 <- computeResiduals(y, X1, theta1)

# --- Model 2: y ~ X4 + X3^2 + X5 ---
X2 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5))
theta2 <- estimate_theta(X2, y)
residuals2 <- computeResiduals(y, X2, theta2)

# --- Model 3: y ~ X3 + X4 + X5^3 ---
X3 <- as.matrix(cbind(1, data$X3, data$X4, data$X5^3))
theta3 <- estimate_theta(X3, y)
residuals3 <- computeResiduals(y, X3, theta3)

# --- Model 4: y ~ X4 + X3^2 + X5^3 ---
X4 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5^3))
theta4 <- estimate_theta(X4, y)
residuals4 <- computeResiduals(y, X4, theta4)

# --- Model 5: y ~ X4 + X1^2 + X3^2 ---
X5 <- as.matrix(cbind(1, data$X4, data$X1^2, data$X3^2))
theta5 <- estimate_theta(X5, y)
residuals5 <- computeResiduals(y, X5, theta5)

# --- Plot Q-Q Plots for Residuals ---
# Function to create Q-Q plot
qq_plot <- function(residuals, title) {
  ggplot(data.frame(residuals), aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "blue") +
    labs(title = title, x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
}

# Generate Q-Q plots for each model
plot1 <- qq_plot(residuals1, "Q-Q Plot for Model 1 Residuals")
plot2 <- qq_plot(residuals2, "Q-Q Plot for Model 2 Residuals")
plot3 <- qq_plot(residuals3, "Q-Q Plot for Model 3 Residuals")
plot4 <- qq_plot(residuals4, "Q-Q Plot for Model 4 Residuals")
plot5 <- qq_plot(residuals5, "Q-Q Plot for Model 5 Residuals")

# Save Q-Q plots
ggsave("E:\\Rpractice\\regression\\qq_plot_model1.png", plot1, width = 6, height = 4)
ggsave("E:\\Rpractice\\regression\\qq_plot_model2.png", plot2, width = 6, height = 4)
ggsave("E:\\Rpractice\\regression\\qq_plot_model3.png", plot3, width = 6, height = 4)
ggsave("E:\\Rpractice\\regression\\qq_plot_model4.png", plot4, width = 6, height = 4)
ggsave("E:\\Rpractice\\regression\\qq_plot_model5.png", plot5, width = 6, height = 4)

# Display all plots
print(plot1)
print(plot2)
print(plot3)
print(plot4)
print(plot5)
