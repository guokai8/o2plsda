################################################################################
# Example: How to plot tuning results
################################################################################

library(o2plsda)
set.seed(123)

# Generate data
X <- scale(matrix(rnorm(5000), 50, 100))
group <- rep(c("Ctrl", "Treat"), each = 25)

# Run tuning
test_values <- c(10, 20, 30, 40, 50)
cat("Running tune_sparse_keepX...\n")
tune_result <- tune_sparse_keepX(X, group,
                                 nc = 2,
                                 test.keepX = test_values,
                                 validation = "Mfold",
                                 folds = 3,
                                 nrepeat = 2)

cat("\nTuning Results:\n")
cat("  Optimal keepX:", tune_result$choice.keepX, "\n")
cat("  Best error rate:", round(tune_result$error.rate, 3), "\n")
cat("\nError rates for all tested values:\n")
print(tune_result$error.rate.all)

# Load plotting functions
source("plot_tune_results.R")

# Create visualizations
cat("\nCreating plots...\n")
p1 <- plot_tune_keepX_results(tune_result, test_values)
p2 <- plot_tune_keepX_barplot(tune_result, test_values)

# Display plots
cat("\nDisplaying line plot...\n")
print(p1)

cat("\nDisplaying barplot...\n")
print(p2)

# Save plots
ggsave("tune_results_line.pdf", p1, width = 8, height = 6)
ggsave("tune_results_bar.pdf", p2, width = 8, height = 6)
cat("\nPlots saved as:\n")
cat("  - tune_results_line.pdf\n")
cat("  - tune_results_bar.pdf\n")

cat("\nDone!\n")
