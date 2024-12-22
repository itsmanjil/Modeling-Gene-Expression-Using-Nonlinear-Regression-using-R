# Load necessary libraries
library(ggplot2)
library(MASS)

# Load the dataset
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine datasets
data <- cbind(time_data, gene_data)
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Assuming X2 is Y

# Split Data into Training (70%) and Testing (30%) ---------------------------
set.seed(123)  # For reproducibility
train_index <- sample(1:nrow(data), 0.7 * nrow(data))  # Random indices for training
train_data <- data[train_index, ]  # Training set
test_data <- data[-train_index, ]  # Testing set

# Select the Best Model: Model 4 (as an example)
# Model 4: y ~ X4 + X3^2 + X5^3
y_train <- as.matrix(train_data$Y)  # Training output
X_train <- as.matrix(cbind(1, train_data$X4, train_data$X3^2, train_data$X5^3))  # Design matrix

# Estimate Parameters using Training Data ------------------------------------
theta <- solve(t(X_train) %*% X_train) %*% t(X_train) %*% y_train
cat("Estimated Parameters (Model 4):", theta, "\n")

# Make Predictions on Testing Data ------------------------------------------
y_test <- as.matrix(test_data$Y)  # Actual output for testing
X_test <- as.matrix(cbind(1, test_data$X4, test_data$X3^2, test_data$X5^3))  # Design matrix for testing
y_pred <- X_test %*% theta  # Predicted output

# Compute Residual Standard Error -------------------------------------------
residuals <- y_train - X_train %*% theta
sigma2 <- sum(residuals^2) / (nrow(X_train) - ncol(X_train))  # Variance estimate
se_pred <- sqrt(sigma2 * diag(X_test %*% solve(t(X_train) %*% X_train) %*% t(X_test)))  # Standard errors

# Compute 95% Confidence Intervals -----------------------------------------
lower_bound <- y_pred - 1.96 * se_pred  # Lower bound
upper_bound <- y_pred + 1.96 * se_pred  # Upper bound

# Combine Predictions and Confidence Intervals
results <- data.frame(
  Time = test_data$Time,
  Actual = y_test,
  Predicted = y_pred,
  Lower = lower_bound,
  Upper = upper_bound
)

# Plot the Predictions with 95% Confidence Intervals -----------------------
plot <- ggplot(results, aes(x = Time)) +
  geom_point(aes(y = Actual), color = "red", size = 2, shape = 16) +  # Testing data
  geom_line(aes(y = Predicted), color = "blue", size = 1) +          # Predicted values
  geom_ribbon(aes(ymin = Lower, ymax = Upper), alpha = 0.3, fill = "blue") +  # Confidence interval
  labs(
    title = "Model Predictions with 95% Confidence Intervals",
    x = "Time",
    y = "Output (Y)"
  ) +
  theme_minimal()

# Save and Display the Plot
ggsave("E:\\Rpractice\\regression\\model_predictions_confidence_intervals.png", plot, width = 10, height = 6)
print(plot)
