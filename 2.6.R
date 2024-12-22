# Load necessary libraries
library(MASS)
library(ggplot2)
library(gridExtra)

# Load the dataset
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine datasets
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Assuming X2 is Y

# Functions ---------------------------------------------------
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

# Function to compute AIC and BIC
computeAIC_BIC <- function(log_likelihood, k, n) {
  AIC <- 2 * k - 2 * log_likelihood
  BIC <- k * log(n) - 2 * log_likelihood
  return(list(AIC = AIC, BIC = BIC))
}

# Function to create Q-Q plot
qq_plot <- function(residuals, title) {
  ggplot(data.frame(residuals), aes(sample = residuals)) +
    stat_qq() +
    stat_qq_line(color = "blue") +
    labs(title = title, x = "Theoretical Quantiles", y = "Sample Quantiles") +
    theme_minimal()
}

# Initialize variables
y <- as.matrix(data$Y)
n <- nrow(data)

# Model Evaluation --------------------------------------------
results <- list()

# --- Model 1: y ~ X4 + X3^2 ---
X1 <- as.matrix(cbind(1, data$X4, data$X3^2))
theta1 <- estimate_theta(X1, y)
rss1 <- computeRSS(y, X1, theta1)
logLik1 <- computeLogLikelihood(n, rss1)
metrics1 <- computeAIC_BIC(logLik1, ncol(X1), n)
residuals1 <- y - X1 %*% theta1
results$Model1 <- list(AIC = metrics1$AIC, BIC = metrics1$BIC, residuals = residuals1)

# --- Model 2: y ~ X4 + X3^2 + X5 ---
X2 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5))
theta2 <- estimate_theta(X2, y)
rss2 <- computeRSS(y, X2, theta2)
logLik2 <- computeLogLikelihood(n, rss2)
metrics2 <- computeAIC_BIC(logLik2, ncol(X2), n)
residuals2 <- y - X2 %*% theta2
results$Model2 <- list(AIC = metrics2$AIC, BIC = metrics2$BIC, residuals = residuals2)

# --- Model 3: y ~ X3 + X4 + X5^3 ---
X3 <- as.matrix(cbind(1, data$X3, data$X4, data$X5^3))
theta3 <- estimate_theta(X3, y)
rss3 <- computeRSS(y, X3, theta3)
logLik3 <- computeLogLikelihood(n, rss3)
metrics3 <- computeAIC_BIC(logLik3, ncol(X3), n)
residuals3 <- y - X3 %*% theta3
results$Model3 <- list(AIC = metrics3$AIC, BIC = metrics3$BIC, residuals = residuals3)

# --- Model 4: y ~ X4 + X3^2 + X5^3 ---
X4 <- as.matrix(cbind(1, data$X4, data$X3^2, data$X5^3))
theta4 <- estimate_theta(X4, y)
rss4 <- computeRSS(y, X4, theta4)
logLik4 <- computeLogLikelihood(n, rss4)
metrics4 <- computeAIC_BIC(logLik4, ncol(X4), n)
residuals4 <- y - X4 %*% theta4
results$Model4 <- list(AIC = metrics4$AIC, BIC = metrics4$BIC, residuals = residuals4)

# --- Model 5: y ~ X4 + X1^2 + X3^2 ---
X5 <- as.matrix(cbind(1, data$X4, data$X1^2, data$X3^2))
theta5 <- estimate_theta(X5, y)
rss5 <- computeRSS(y, X5, theta5)
logLik5 <- computeLogLikelihood(n, rss5)
metrics5 <- computeAIC_BIC(logLik5, ncol(X5), n)
residuals5 <- y - X5 %*% theta5
results$Model5 <- list(AIC = metrics5$AIC, BIC = metrics5$BIC, residuals = residuals5)

# Compare Results ---------------------------------------------
cat("AIC and BIC Comparison:\n")
for (model in names(results)) {
  cat(model, "AIC:", results[[model]]$AIC, "BIC:", results[[model]]$BIC, "\n")
}

# Generate Q-Q plots for all models
plot1 <- qq_plot(results$Model1$residuals, "Q-Q Plot for Model 1 Residuals")
plot2 <- qq_plot(results$Model2$residuals, "Q-Q Plot for Model 2 Residuals")
plot3 <- qq_plot(results$Model3$residuals, "Q-Q Plot for Model 3 Residuals")
plot4 <- qq_plot(results$Model4$residuals, "Q-Q Plot for Model 4 Residuals")
plot5 <- qq_plot(results$Model5$residuals, "Q-Q Plot for Model 5 Residuals")

# Display Q-Q plots
grid.arrange(plot1, plot2, plot3, plot4, plot5, ncol = 2)
