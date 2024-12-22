# Load necessary libraries
library(ggplot2)
library(reshape2)

# Load the datasets
gene_data <- read.csv("E:\\Rpractice\\regression\\data_03cf4b0d-8941-4323-998c-a7ff1a83d0f2_1733812247131.csv")
time_data <- read.csv("E:\\Rpractice\\regression\\time_1673241270748.csv")

# Combine the time and gene expression data
data <- cbind(time_data, gene_data)

# Rename the columns for clarity
colnames(data) <- c("Time", "X1", "Y", "X3", "X4", "X5")  # Assuming X2 is Y output

# --- 1. Time Series Plots for X Inputs (X1, X3, X4, X5) ---
# Melt the data to long format for plotting
data_long <- melt(data, id.vars = "Time", measure.vars = c("X1", "X3", "X4", "X5"),
                  variable.name = "Gene", value.name = "Expression")

# Generate time series plot for X inputs
plot_x <- ggplot(data_long, aes(x = Time, y = Expression, color = Gene)) +
  geom_line(size = 0.8) +
  facet_wrap(~Gene, scales = "free_y", ncol = 1) +
  labs(title = "Time Series Analysis of X Inputs",
       x = "Time", y = "Expression Levels") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue", color = "white"),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none"
  )

ggsave("E:\\Rpractice\\regression\\time_series_X_inputs.png", plot_x, width = 10, height = 8)

# --- 2. Time Series Plot for Y Output ---
plot_y <- ggplot(data, aes(x = Time, y = Y)) +
  geom_line(color = "steelblue", size = 1) +
  labs(title = "Time Series Analysis of Y with Time",
       x = "Time", y = "Output Signal Y") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue", color = "white"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Save Y output plot
ggsave("E:\\Rpractice\\regression\\time_series_Y_output.png", plot_y, width = 10, height = 6)
# Generate density distribution plot
plot_density <- ggplot(data_long, aes(x = Expression, fill = Gene, color = Gene)) +
  geom_histogram(aes(y = ..density..), bins = 30, alpha = 0.4, position = "identity") +
  geom_density(size = 1) +
  labs(title = "Distribution of X Inputs",
       x = "X Signal", y = "Density") +
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "lightblue", color = "white"),
    legend.title = element_text(size = 12, face = "bold"),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

ggsave("E:\\Rpractice\\regression\\density_distribution_X_inputs.png", plot_density, width = 10, height = 6)
# Generate individual density plots (excluding X2)
plot_X1 <- ggplot(data, aes(x = X1)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "white") +
  geom_density(color = "black", size = 1) +
  labs(title = "Distribution of X1", x = "X1 Signal", y = "Density") +
  theme_minimal()

plot_X3 <- ggplot(data, aes(x = X3)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "white") +
  geom_density(color = "black", size = 1) +
  labs(title = "Distribution of X3", x = "X3 Signal", y = "Density") +
  theme_minimal()

plot_X4 <- ggplot(data, aes(x = X4)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "white") +
  geom_density(color = "black", size = 1) +
  labs(title = "Distribution of X4", x = "X4 Signal", y = "Density") +
  theme_minimal()

plot_X5 <- ggplot(data, aes(x = X5)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "white") +
  geom_density(color = "black", size = 1) +
  labs(title = "Distribution of X5", x = "X5 Signal", y = "Density") +
  theme_minimal()

# Arrange the plots in a 2x2 grid
grid_plot <- grid.arrange(plot_X1, plot_X3, plot_X4, plot_X5, nrow = 2, ncol = 2)

# Save the plot as an image file
ggsave("E:\\Rpractice\\regression\\individual_X_density_plots.png", grid_plot, width = 12, height = 8)

# Generate the distribution plot for Y
plot_Y <- ggplot(data, aes(x = Y)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", color = "white") +
  geom_density(color = "black", size = 1) +
  labs(title = "Distribution of Y",
       x = "Y Signal", y = "Density") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  )

# Save the plot
ggsave("E:\\Rpractice\\regression\\distribution_Y_output.png", plot_Y, width = 10, height = 6)

# Melt the data to long format for scatter plotting
data_long <- melt(data, id.vars = "Y", measure.vars = c("X1", "X3", "X4", "X5"),
                  variable.name = "X_inputs", value.name = "X_values")

# Generate scatter plot for correlation
correlation_plot <- ggplot(data_long, aes(x = X_values, y = Y, color = X_inputs)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  labs(title = "Correlation of X Inputs with Y Output",
       x = "X Inputs", y = "Y") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.title = element_blank()
  )

# Save the plot
ggsave("E:\\Rpractice\\regression\\correlation_X_Y_plot.png", correlation_plot, width = 10, height = 6)

# Generate individual scatter plots
plot_X1 <- ggplot(data, aes(x = X1, y = Y)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  labs(title = "Correlation of X1 input to Y output",
       x = "X1 input", y = "Output Signal Y") +
  theme_minimal()

plot_X3 <- ggplot(data, aes(x = X3, y = Y)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  labs(title = "Correlation of X3 input to Y output",
       x = "X3 input", y = "Output Signal Y") +
  theme_minimal()

plot_X4 <- ggplot(data, aes(x = X4, y = Y)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  labs(title = "Correlation of X4 input to Y output",
       x = "X4 input", y = "Output Signal Y") +
  theme_minimal()

plot_X5 <- ggplot(data, aes(x = X5, y = Y)) +
  geom_point(color = "blue") +
  geom_hline(yintercept = 0, color = "black") +
  geom_vline(xintercept = 0, color = "black") +
  labs(title = "Correlation of X5 input to Y output",
       x = "X5 input", y = "Output Signal Y") +
  theme_minimal()

# Arrange the plots in a 2x2 grid
grid_plot <- grid.arrange(plot_X1, plot_X3, plot_X4, plot_X5, nrow = 2, ncol = 2)

# Save the plot as an image file
ggsave("E:\\Rpractice\\regression\\correlation_X_individual_Y.png", grid_plot, width = 12, height = 8)

# Function to generate Model 1
generateModel1 <- function(df) {
  set.seed(100)
  ones <- matrix(1, length(df$X1), 1)
  noise <- rnorm(length(df$y), mean = 0, sd = 0.1)
  theta_bias <- runif(1, -1, 1) * ones
  return(cbind(df$X4, df$X1^2, df$X1^3, df$X2^4, df$X1^4, theta_bias, noise))
}

# Function to generate Model 2
generateModel2 <- function(df) {
  set.seed(100)
  ones <- matrix(1, length(df$X1), 1)
  noise <- rnorm(length(df$y), mean = 0, sd = 0.1)
  theta_bias <- runif(1, -1, 1) * ones
  return(cbind(df$X4, df$X1^3, df$X3^4, theta_bias, noise))
}

# Function to generate Model 3
generateModel3 <- function(df) {
  set.seed(100)
  ones <- matrix(1, length(df$X1), 1)
  noise <- rnorm(length(df$y), mean = 0, sd = 0.1)
  theta_bias <- runif(1, -1, 1) * ones
  return(cbind(df$X3^3, df$X3^4, theta_bias, noise))
}

# Function to generate Model 4
generateModel4 <- function(df) {
  set.seed(100)
  ones <- matrix(1, length(df$X1), 1)
  noise <- rnorm(length(df$y), mean = 0, sd = 0.1)
  theta_bias <- runif(1, -1, 1) * ones
  return(cbind(df$X2, df$X1^3, df$X3^4, theta_bias, noise))
}

# Function to generate Model 5
generateModel5 <- function(df) {
  set.seed(100)
  ones <- matrix(1, length(df$X1), 1)
  noise <- rnorm(length(df$y), mean = 0, sd = 0.1)
  theta_bias <- runif(1, -1, 1) * ones
  return(cbind(df$X4, df$X1^2, df$X1^3, df$X3^4, theta_bias, noise))
}


# Model 1: y ~ X4 + X3^2
model1 <- lm(Y ~ X4 + I(X3^2), data = data)
cat("Model 1 Parameters:\n")
print(coef(model1))

# Model 2: y ~ X4 + X3^2 + X5
model2 <- lm(Y ~ X4 + I(X3^2) + X5, data = data)
cat("\nModel 2 Parameters:\n")
print(coef(model2))

# Model 3: y ~ X3 + X4 + X5^3
model3 <- lm(Y ~ X3 + X4 + I(X5^3), data = data)
cat("\nModel 3 Parameters:\n")
print(coef(model3))

# Model 4: y ~ X4 + X3^2 + X5^3
model4 <- lm(Y ~ X4 + I(X3^2) + I(X5^3), data = data)
cat("\nModel 4 Parameters:\n")
print(coef(model4))

# Model 5: y ~ X4 + X1^2 + X3^2
model5 <- lm(Y ~ X4 + I(X1^2) + I(X3^2), data = data)
cat("\nModel 5 Parameters:\n")
print(coef(model5))



# Display success message
print("Individual correlation plots for X inputs with Y output saved successfully!")
print("Correlation plot for X inputs and Y output saved successfully to 'E:\\Rpractice\\regression\\correlation_X_Y_plot.png'.")
print("Distribution plot for Y saved successfully to 'E:\\Rpractice\\regression\\distribution_Y_output.png'.")
print("Individual density plots for X1, X3, X4, and X5 saved successfully to 'E:\\Rpractice\\regression\\individual_X_density_plots.png'.")
print("Density plot for X inputs has been saved to 'E:\\Rpractice\\regression\\density_distribution_X_inputs.png'.")
print("X inputs time series plot has been saved to 'E:\\Rpractice\\regression\\time_series_X_inputs.png'.")
print("Y output time series plot has been saved to 'E:\\Rpractice\\regression\\time_series_Y_output.png'.")
