# Comprehensive Tutorial Guide for o2plsda Package

## Overview

This tutorial provides a complete walkthrough of **all major functions** in the o2plsda package. It's designed to help you understand and test every feature available in the package.

## Tutorial Structure

The tutorial is organized into 15 parts, each focusing on different aspects of the package:

### Part 1: Data Preparation
- Generate example omics data matrices (X and Y)
- Create sample metadata and group labels
- Data scaling and normalization

### Part 2: Standard O2PLS Analysis
- **o2cv()**: Cross-validation for parameter selection
- **o2pls()**: Fit O2PLS model
- **summary()** and **print()**: Examine model results
- **loadings()**: Extract loading matrices (Xjoint, Yjoint, Xorth, Yorth)
- **scores()**: Extract score matrices

### Part 3: O2PLS Plotting
- Score plots for joint and orthogonal components
- Loading plots with various options
- Rotation and repel options

### Part 4: O2PLS-DA (Discriminant Analysis)
- **oplsda()**: Perform OPLS-DA on O2PLS results
- **vip()**: Calculate Variable Importance in Projection
- Visualization of discriminant analysis results

### Part 5: Standard PLS-DA
- **plsda()**: Direct PLS-DA analysis
- Extract components and VIP scores
- Score, VIP, and loading plots

### Part 6: Sparse O2PLS
- **sparse_o2pls()**: Feature selection with sparsity
- **selected_vars()**: Get indices of selected variables
- **selected_var_names()**: Get names of selected variables
- **sparsity_info()**: Get sparsity statistics
- **predict()**: Make predictions on new data
- Multiple plot types (score, loading, biplot, diagnostic, contribution)

### Part 7: Sparse PLS-DA
- **sparse_plsda()**: Sparse discriminant analysis
- Variable selection with cross-validation
- Classification and predictions
- Comprehensive visualization options

### Part 8: Tuning Sparse Models
- **tune_sparse_keepX()**: Optimize number of variables to keep
- **tune_sparse_o2pls()**: Optimize all sparse O2PLS parameters
- Visualization of tuning results (heatmap, line plots)
- **Note:** Uses custom plotting functions from `plot_tune_results.R` for tune_sparse_keepX results

### Part 9: Cross-Validation for Sparse O2PLS
- **o2cv_sparse()**: Cross-validation with sparsity constraints
- Parameter optimization

### Part 10: Stability Selection
- **stability_selection()**: Bootstrap-based feature selection
- Multiple visualization types:
  - Heatmap of selection frequencies
  - Barplot of stable features
  - Network plot of feature relationships
  - Summary plots

### Part 11: Model Comparison
- **compare_models()**: Compare multiple fitted models
- Performance metrics comparison

### Part 12: Advanced Prediction Functions
- **calculate_mahalanobis_distances()**: Distance-based classification
- **calculate_prediction_probabilities()**: Probability estimates

### Part 13: Optimal Parameter Finding
- **find_optimal_sparse_params()**: Automated parameter optimization
- Grid search over parameter space

### Part 14: All Sparse O2PLS Plot Types
Demonstrates all available plot types for Sparse O2PLS:
- Score plots (Xjoint, Yjoint)
- Loading plots (Xjoint, Yjoint)
- Biplots (X, Y)
- Diagnostic plots (X, Y)
- Contribution plots (X, Y)

### Part 15: All Sparse PLS-DA Plot Types
Demonstrates all available plot types for Sparse PLS-DA:
- Score plot
- Loading plot
- VIP plot
- Contribution plot
- Selection plot
- Classification plot
- Diagnostic plot

## How to Run the Tutorial

### Option 1: Run the entire script
```r
source("comprehensive_tutorial.R")
```

### Option 2: Run section by section
Open the file in RStudio and run each section individually using Ctrl+Enter (Windows/Linux) or Cmd+Enter (Mac).

### Option 3: Run specific parts
```r
# Load the library
library(o2plsda)

# Copy and paste specific sections you want to test
# For example, to test only Sparse O2PLS:
# Run PART 1 (Data Preparation)
# Then run PART 6 (Sparse O2PLS)
```

## Expected Output

The tutorial will:
1. Print progress messages for each step
2. Display summaries of fitted models
3. Show top selected variables and VIP scores
4. Create 41+ different plots (stored in R environment as p1, p2, ..., p41)
5. Generate predictions and probability estimates

## Viewing Plots

All plots are stored as ggplot objects. To view them:

```r
# View individual plots
print(p1)  # First O2PLS score plot
print(p15) # Sparse PLS-DA score plot

# View plots in a grid
library(gridExtra)
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Save plots to files
ggsave("o2pls_score_plot.pdf", p1, width = 8, height = 6)
```

## Functions Tested

The tutorial comprehensively tests **24 main functions**:

### Core Analysis Functions
1. `o2pls()` - O2PLS analysis
2. `oplsda()` - O2PLS discriminant analysis
3. `plsda()` - PLS discriminant analysis
4. `sparse_o2pls()` - Sparse O2PLS
5. `sparse_plsda()` - Sparse PLS-DA

### Cross-Validation & Tuning
6. `o2cv()` - O2PLS cross-validation
7. `o2cv_sparse()` - Sparse O2PLS cross-validation
8. `tune_sparse_keepX()` - Tune sparsity parameter
9. `tune_sparse_o2pls()` - Tune all parameters
10. `find_optimal_sparse_params()` - Automated optimization

### Feature Selection
11. `stability_selection()` - Stability-based selection
12. `selected_vars()` - Get selected variable indices
13. `selected_var_names()` - Get selected variable names
14. `sparsity_info()` - Get sparsity statistics

### Extraction Functions
15. `loadings()` - Extract loadings
16. `scores()` - Extract scores
17. `vip()` - Variable importance

### Prediction & Classification
18. `predict()` - Make predictions
19. `calculate_mahalanobis_distances()` - Distance calculation
20. `calculate_prediction_probabilities()` - Probability calculation

### Utilities
21. `compare_models()` - Model comparison
22. `summary()` - Model summaries
23. `print()` - Print model info
24. `plot()` - Comprehensive visualization

## Plot Types Available

### O2PLS Plots
- `type = "score"` with `var = "Xjoint"/"Yjoint"/"Xorth"/"Yorth"`
- `type = "loading"` with `var = "Xjoint"/"Yjoint"/"Xorth"/"Yorth"`

### O2PLS-DA Plots
- `type = "score"`
- `type = "vip"`
- `type = "loading"`

### PLS-DA Plots
- `type = "score"`
- `type = "vip"`
- `type = "loading"`

### Sparse O2PLS Plots
- `type = "score"`
- `type = "loading"`
- `type = "biplot"`
- `type = "diagnostic"`
- `type = "contribution"`

### Sparse PLS-DA Plots
- `type = "score"`
- `type = "loading"`
- `type = "vip"`
- `type = "contribution"`
- `type = "selection"`
- `type = "classification"`
- `type = "diagnostic"`

### Tuning Result Plots
- `type = "heatmap"`
- `type = "line"`

### Stability Selection Plots
- `type = "heatmap"`
- `type = "barplot"`
- `type = "network"`
- `type = "summary"`

## Customization Options

### Plot Customization
```r
# Add custom colors, themes, labels, etc.
p <- plot(o2pls_fit, type = "score", var = "Xjoint", group = group)
p <- p + ggplot2::theme_bw() + ggplot2::ggtitle("My Custom Title")
```

### Model Parameters
```r
# Adjust model parameters based on your data
o2pls_fit <- o2pls(X, Y,
                   nc = 3,    # More joint components
                   nx = 1,    # Fewer X-specific components
                   ny = 1)    # Fewer Y-specific components
```

### Sparsity Levels
```r
# Control feature selection
sparse_fit <- sparse_o2pls(X, Y, nc = 1, nx = 0, ny = 0,
                           keepX = 50,   # Keep more X features
                           keepY = 10)   # Keep fewer Y features
```

## Tips for Your Own Data

1. **Data Scaling**: Always scale your data before analysis
   ```r
   X <- scale(X, center = TRUE, scale = TRUE)
   ```

2. **Parameter Selection**: Use cross-validation to find optimal parameters
   ```r
   cv <- o2cv(X, Y, 1:5, 0:3, 0:3, nr_folds = 10)
   ```

3. **Sparsity**: Start with a reasonable percentage of features (e.g., 20-30%)
   ```r
   keepX <- round(ncol(X) * 0.3)  # Keep 30% of features
   ```

4. **Sample Size**: Ensure you have enough samples for cross-validation
   - Minimum: 20-30 samples per group
   - Recommended: 50+ samples per group

5. **Feature Size**: The package handles high-dimensional data well
   - Works with p >> n (more features than samples)
   - Use sparse methods for very high dimensions

## Troubleshooting

### Common Issues

1. **"Matrix not full rank" error**
   - Solution: Reduce number of components or increase sample size

2. **Slow cross-validation**
   - Solution: Use parallel processing with `nr_cores` parameter
   ```r
   cv <- o2cv(X, Y, 1:3, 0:2, 0:2, nr_cores = 4)
   ```

3. **Memory issues with large data**
   - Solution: Use sparse methods and limit keepX/keepY values

4. **Plot labels overlapping**
   - Solution: Use `repel = TRUE` for automatic label adjustment
   ```r
   plot(fit, type = "loading", repel = TRUE)
   ```

5. **"'x' is a list, but does not have components 'x' and 'y'" when plotting tune_result**
   - Issue: `tune_sparse_keepX()` returns a simple list, not a TuneResult object
   - Solution: Use custom plotting functions from `plot_tune_results.R`
   ```r
   # Instead of: plot(tune_result, type = "heatmap")
   # Use:
   source("plot_tune_results.R")
   p <- plot_tune_keepX_results(tune_result, test.keepX)
   ```
   - See `PLOTTING_TUNING_RESULTS.md` for full details

## Performance Optimization

For large datasets:

```r
# Use parallel processing
cv <- o2cv(X, Y, 1:3, 0:2, 0:2,
           nr_cores = parallel::detectCores() - 1)

# Reduce cross-validation iterations
tune <- tune_sparse_keepX(X, group, nc = 2,
                          folds = 3,      # Instead of 10
                          nrepeat = 2)    # Instead of 10

# Use stability selection with fewer bootstraps
stab <- stability_selection(X, group, nc = 2,
                            n_bootstrap = 20)  # Instead of 100
```

## Saving Results

```r
# Save all results
save(o2pls_fit, oplsda_result, sparse_o2pls_fit,
     file = "my_analysis_results.RData")

# Save plots
pdf("all_plots.pdf", width = 10, height = 8)
print(p1)
print(p2)
# ... etc
dev.off()

# Export selected variables
write.csv(data.frame(
  variable = selected_var_names(sparse_plsda_fit),
  vip = vip_scores[selected_vars(sparse_plsda_fit)]
), file = "selected_variables.csv")
```

## Next Steps

After running this tutorial:

1. **Analyze your own data** using the patterns shown
2. **Customize visualizations** to match your publication needs
3. **Export results** for further analysis or reporting
4. **Consult the package vignette** for additional details
5. **Check the help files** for individual functions: `?o2pls`, `?sparse_plsda`, etc.

## References

- Trygg J, Wold S (2003). "O2-PLS, a two-block (X-Y) latent variable regression (LVR) method with an integral OSC filter." *Journal of Chemometrics*, 17(1), 53-64.
- Löfstedt T, Trygg J (2011). "OnPLS—a novel multiblock method for the modelling of predictive and orthogonal variation." *Journal of Chemometrics*, 25(8), 441-455.

## Support

For questions, issues, or feature requests:
- GitHub: https://github.com/guokai8/o2plsda/issues
- Email: guokai8@gmail.com

---

**Last Updated**: December 2025
**Package Version**: 0.0.28
