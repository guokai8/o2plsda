################################################################################
# Workaround: Visualize tune_sparse_keepX results
# The tune_result from tune_sparse_keepX is a simple list, not a TuneResult object
################################################################################

library(ggplot2)

#' Plot tuning results from tune_sparse_keepX
#'
#' @param tune_result Result from tune_sparse_keepX()
#' @param test.keepX Vector of keepX values tested
#' @return ggplot object
plot_tune_keepX_results <- function(tune_result, test.keepX) {

  # Extract error rates
  error_rates <- tune_result$error.rate.all

  # Determine if it's matrix (multiple components) or vector (single component)
  if(is.matrix(error_rates)) {
    # Multiple components case
    nc <- ncol(error_rates)

    # Create data frame
    df <- data.frame()
    for(comp in 1:nc) {
      for(i in 1:length(test.keepX)) {
        df <- rbind(df, data.frame(
          keepX = test.keepX[i],
          component = paste0("Comp_", comp),
          error_rate = error_rates[i, comp]
        ))
      }
    }

    # Plot as line plot with facets
    p <- ggplot(df, aes(x = keepX, y = error_rate, color = component, group = component)) +
      geom_line(linewidth = 1) +
      geom_point(size = 3) +
      labs(x = "Number of Variables (keepX)",
           y = "Error Rate",
           title = "Tuning Results: Error Rate vs keepX",
           subtitle = paste("Optimal keepX:", paste(tune_result$choice.keepX, collapse = ", "))) +
      theme_minimal() +
      theme(legend.position = "bottom")

  } else {
    # Single component case
    df <- data.frame(
      keepX = test.keepX,
      error_rate = as.vector(error_rates)
    )

    # Find optimal
    optimal_idx <- which(test.keepX == tune_result$choice.keepX)

    # Plot as simple line plot
    p <- ggplot(df, aes(x = keepX, y = error_rate)) +
      geom_line(linewidth = 1, color = "steelblue") +
      geom_point(size = 3, color = "steelblue") +
      geom_point(data = df[optimal_idx, ], aes(x = keepX, y = error_rate),
                 size = 5, color = "red", shape = 18) +
      geom_vline(xintercept = tune_result$choice.keepX,
                 linetype = "dashed", color = "red", alpha = 0.5) +
      annotate("text", x = tune_result$choice.keepX,
               y = max(df$error_rate) * 0.95,
               label = paste("Optimal =", tune_result$choice.keepX),
               hjust = -0.1, color = "red") +
      labs(x = "Number of Variables (keepX)",
           y = "Error Rate",
           title = "Tuning Results: Error Rate vs keepX",
           subtitle = paste("Optimal keepX:", tune_result$choice.keepX,
                          "| Error Rate:", round(tune_result$error.rate, 3))) +
      theme_minimal()
  }

  return(p)
}

#' Alternative: Simple barplot
plot_tune_keepX_barplot <- function(tune_result, test.keepX) {

  error_rates <- tune_result$error.rate.all

  if(is.matrix(error_rates)) {
    # Take mean across components
    mean_errors <- rowMeans(error_rates, na.rm = TRUE)
  } else {
    mean_errors <- as.vector(error_rates)
  }

  df <- data.frame(
    keepX = factor(test.keepX),
    error_rate = mean_errors,
    is_optimal = test.keepX == tune_result$choice.keepX[1]
  )

  p <- ggplot(df, aes(x = keepX, y = error_rate, fill = is_optimal)) +
    geom_col() +
    scale_fill_manual(values = c("TRUE" = "darkred", "FALSE" = "steelblue"),
                      guide = "none") +
    geom_text(aes(label = round(error_rate, 3)), vjust = -0.5, size = 3) +
    labs(x = "Number of Variables (keepX)",
         y = "Mean Error Rate",
         title = "Tuning Results",
         subtitle = "Red bar indicates optimal choice") +
    theme_minimal() +
    theme(panel.grid.major.x = element_blank())

  return(p)
}

################################################################################
# Example usage:
################################################################################

# After running tune_sparse_keepX:
# library(o2plsda)
# set.seed(123)
# X <- scale(matrix(rnorm(5000), 50, 100))
# group <- rep(c("Ctrl", "Treat"), each = 25)
#
# test_values <- c(10, 20, 30, 40, 50)
# tune_result <- tune_sparse_keepX(X, group, nc = 2,
#                                  test.keepX = test_values,
#                                  validation = "Mfold",
#                                  folds = 3, nrepeat = 2)
#
# # Now plot:
# source("plot_tune_results.R")
# p1 <- plot_tune_keepX_results(tune_result, test_values)
# print(p1)
#
# p2 <- plot_tune_keepX_barplot(tune_result, test_values)
# print(p2)

cat("Functions loaded:\n")
cat("  - plot_tune_keepX_results(tune_result, test.keepX)\n")
cat("  - plot_tune_keepX_barplot(tune_result, test.keepX)\n")
cat("\nUse these to visualize tune_sparse_keepX results.\n")
