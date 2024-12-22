# Install necessary libraries
if (!require(ggplot2)) install.packages("ggplot2", dependencies = TRUE)

# Plot distribution of training data
training_plot <- ggplot(data = data.frame(y = training_set$Y), aes(x = y)) +
  geom_density(color = "darkgreen", fill = "lightgreen", alpha = 0.4) +
  geom_vline(xintercept = mean(training_set$Y), color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = quantile(training_set$Y, probs = c(0.025, 0.975)), color = "darkgreen", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of training data", x = "y", y = "Density") +
  theme_minimal()

# Save the training plot
ggsave("training_distribution.png", plot = training_plot, width = 8, height = 6, dpi = 300)

# Plot distribution of testing data
testing_plot <- ggplot(data = data.frame(y = testing_set$Y), aes(x = y)) +
  geom_density(color = "blue", fill = "lightblue", alpha = 0.4) +
  geom_vline(xintercept = mean(testing_set$Y), color = "black", linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = quantile(testing_set$Y, probs = c(0.025, 0.975)), color = "blue", linetype = "dashed", linewidth = 1) +
  labs(title = "Distribution of testing data", x = "y", y = "Density") +
  theme_minimal()

# Save the testing plot
ggsave("testing_distribution.png", plot = testing_plot, width = 8, height = 6, dpi = 300)

# Print confirmation
cat("The training distribution plot has been saved as: training_distribution.png\n")
cat("The testing distribution plot has been saved as: testing_distribution.png\n")

# Print plots to the console
print(training_plot)
print(testing_plot)
