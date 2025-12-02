# Plotting Tuning Results - Workaround

## The Issue

When you run `tune_sparse_keepX()` and try to plot the results:

```r
tune_result <- tune_sparse_keepX(X, group, nc = 2, test.keepX = c(10, 20, 30))
plot(tune_result, type = "heatmap")  # ERROR!
```

**Error:**
```
Error in xy.coords(x, y, xlabel, ylabel, log) :
  'x' is a list, but does not have components 'x' and 'y'
```

## Why This Happens

- `tune_sparse_keepX()` returns a **simple list** with elements: `choice.keepX`, `error.rate`, `error.rate.all`
- `plot.TuneResult()` expects a **TuneResult S4 object** (used by other tuning functions like `tune_sparse_o2pls()`)
- The default `plot()` method doesn't know how to handle this simple list format

## The Solution

Use the custom plotting functions provided in `plot_tune_results.R`:

### Quick Usage

```r
# After running tune_sparse_keepX
test_values <- c(10, 20, 30, 40, 50)
tune_result <- tune_sparse_keepX(X, group, nc = 2, test.keepX = test_values)

# Load custom plotting functions
source("plot_tune_results.R")

# Create line plot
p1 <- plot_tune_keepX_results(tune_result, test_values)
print(p1)

# Create barplot
p2 <- plot_tune_keepX_barplot(tune_result, test_values)
print(p2)
```

## Two Plotting Options

### Option 1: Line Plot (Recommended)
Shows error rates across all keepX values tested.

```r
p <- plot_tune_keepX_results(tune_result, test.keepX)
```

**Features:**
- Shows error rate trend
- Highlights optimal choice in red
- Handles multiple components (faceted)
- Includes dashed line at optimal value

### Option 2: Barplot
Simple barplot of error rates.

```r
p <- plot_tune_keepX_barplot(tune_result, test.keepX)
```

**Features:**
- Clear comparison of all values
- Optimal choice highlighted in dark red
- Shows error rate values on top of bars
- Good for presentations

## Complete Example

```r
library(o2plsda)
set.seed(123)

# Prepare data
X <- scale(matrix(rnorm(5000), 50, 100))
group <- rep(c("Ctrl", "Treat"), each = 25)

# Define values to test
test_values <- seq(10, 50, 10)  # Test 10, 20, 30, 40, 50

# Run tuning
tune_result <- tune_sparse_keepX(
  X, group,
  nc = 2,
  test.keepX = test_values,
  validation = "Mfold",
  folds = 5,
  nrepeat = 10
)

# View results
cat("Optimal keepX:", tune_result$choice.keepX, "\n")
cat("Error rate:", round(tune_result$error.rate, 3), "\n")

# Load plotting functions
source("plot_tune_results.R")

# Create and save plots
p1 <- plot_tune_keepX_results(tune_result, test_values)
p2 <- plot_tune_keepX_barplot(tune_result, test_values)

# Display
print(p1)
print(p2)

# Save
ggsave("tuning_line.pdf", p1, width = 8, height = 6)
ggsave("tuning_bar.pdf", p2, width = 6, height = 5)
```

## Understanding the Results

### Interpreting tune_result

The `tune_result` object contains:

1. **choice.keepX** - The optimal number of variables to keep
   ```r
   tune_result$choice.keepX
   # [1] 20  (for single component)
   # [1] 20 30  (for two components - one value per component)
   ```

2. **error.rate** - The minimum error rate achieved
   ```r
   tune_result$error.rate
   # [1] 0.234
   ```

3. **error.rate.all** - Error rates for all tested values
   ```r
   tune_result$error.rate.all
   #      [,1]  [,2]   (if nc=2)
   # [1,] 0.45  0.42
   # [2,] 0.34  0.38
   # [3,] 0.28  0.30   <- lowest for component 1
   ```

### Manual Inspection

```r
# View all error rates
print(tune_result$error.rate.all)

# Find which keepX gave the best result
best_idx <- which.min(tune_result$error.rate.all[, 1])
best_keepX <- test_values[best_idx]

# Plot manually with base R
plot(test_values, tune_result$error.rate.all[, 1],
     type = "b", pch = 19,
     xlab = "keepX", ylab = "Error Rate",
     main = "Tuning Results")
abline(v = tune_result$choice.keepX, col = "red", lty = 2)
```

## Alternative: Use tune_sparse_o2pls() Instead

If you need the full plotting capabilities, use `tune_sparse_o2pls()` which returns a proper TuneResult object:

```r
# This DOES work with plot()
tune_o2pls <- tune_sparse_o2pls(
  X, Y,
  nc_range = 1:3,
  nx_range = 0:2,
  ny_range = 0:2,
  keepX_range = c(20, 30, 40),
  keepY_range = c(20, 30, 40),
  nr_folds = 5
)

# This works!
plot(tune_o2pls, type = "heatmap")
plot(tune_o2pls, type = "line")
```

## Quick Reference

| Function | Returns | plot() works? | Use instead |
|----------|---------|---------------|-------------|
| `tune_sparse_keepX()` | Simple list | ❌ No | `plot_tune_keepX_results()` |
| `tune_sparse_o2pls()` | TuneResult S4 | ✅ Yes | Regular `plot()` ✓ |
| `o2cv_sparse()` | List | ❌ No | Manual visualization |

### Parameter Reference

| Function | Folds Parameter | Cores Parameter |
|----------|-----------------|-----------------|
| `o2cv()` | `nr_folds` | `ncores` |
| `o2cv_sparse()` | `nr_folds` | `ncores` |
| `tune_sparse_keepX()` | `folds` | `nrepeat` (no parallel) |
| `tune_sparse_o2pls()` | `folds` | N/A (no parallel) |

## Files Provided

1. **plot_tune_results.R** - Custom plotting functions
2. **example_tune_plot.R** - Complete working example
3. **PLOTTING_TUNING_RESULTS.md** - This guide

## Support

If you encounter issues:
- Check that you're passing the correct `test.keepX` values to the plotting function
- Ensure the tune_result object has the expected structure
- Try the example script: `source("example_tune_plot.R")`

---

**Last Updated:** December 2025
