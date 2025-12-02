# Parameter Corrections for Tutorial

## Issue Found

The comprehensive tutorial had incorrect parameter names for some functions.

## Corrections Made

### 1. tune_sparse_o2pls()

**Before (INCORRECT):**
```r
tune_o2pls <- tune_sparse_o2pls(X, Y,
                                 nc_range = 1:2,
                                 nx_range = 0:1,
                                 ny_range = 0:1,
                                 keepX_range = c(20, 30),
                                 keepY_range = c(20, 30),
                                 nr_folds = 3,      # ❌ Wrong parameter name
                                 nr_cores = 1)      # ❌ Doesn't exist
```

**After (CORRECT):**
```r
tune_o2pls <- tune_sparse_o2pls(X, Y,
                                 nc_range = 1:2,
                                 nx_range = 0:1,
                                 ny_range = 0:1,
                                 keepX_range = c(20, 30),
                                 keepY_range = c(20, 30),
                                 folds = 3,         # ✅ Correct
                                 measure = "RMSE")  # ✅ Correct
```

**Note:** This function does NOT have parallel processing support.

---

### 2. o2cv_sparse()

**Before (INCORRECT):**
```r
cv_sparse <- o2cv_sparse(X, Y,
                         nc = 1,
                         nx = 1,
                         ny = 1,
                         keepX_range = c(20, 30),
                         keepY_range = c(20, 30),
                         nr_folds = 3,
                         nr_cores = 1)  # ❌ Wrong parameter name
```

**After (CORRECT):**
```r
cv_sparse <- o2cv_sparse(X, Y,
                         nc = 1,
                         nx = 1,
                         ny = 1,
                         keepX_range = c(20, 30),
                         keepY_range = c(20, 30),
                         nr_folds = 3,
                         ncores = 1)    # ✅ Correct (ncores not nr_cores)
```

---

### 3. stability_selection()

**Before (INCORRECT):**
```r
stability_result <- stability_selection(X, group,
                                         nc = 2,
                                         keepX_range = seq(10, 40, 10),
                                         n_bootstrap = 20,
                                         folds = 3,      # ❌ Doesn't exist
                                         threshold = 0.6)
```

**After (CORRECT):**
```r
stability_result <- stability_selection(X, group,
                                         nc = 2,
                                         keepX_range = seq(10, 40, 10),
                                         n_bootstrap = 20,
                                         threshold = 0.6,
                                         parallel = FALSE)  # ✅ Correct
```

**Note:** This function does NOT have a `folds` parameter. It uses bootstrap resampling instead.

---

## Complete Parameter Reference

### Cross-Validation Functions

| Function | Folds Parameter | Parallel Parameter | Notes |
|----------|-----------------|-------------------|-------|
| `o2cv()` | `nr_folds` | `ncores` | Standard O2PLS CV |
| `o2cv_sparse()` | `nr_folds` | `ncores` | Sparse O2PLS CV |

### Tuning Functions

| Function | Folds Parameter | Parallel Parameter | Notes |
|----------|-----------------|-------------------|-------|
| `tune_sparse_keepX()` | `folds` | None (uses `nrepeat`) | Returns simple list |
| `tune_sparse_o2pls()` | `folds` | None | Returns TuneResult S4 |

### Selection Functions

| Function | Bootstrap Parameter | Parallel Parameter | Notes |
|----------|---------------------|-------------------|-------|
| `stability_selection()` | `n_bootstrap` | `parallel`, `ncores` | No folds parameter |

---

## Key Differences

### Folds Parameter Names
- **`nr_folds`**: Used by `o2cv()` and `o2cv_sparse()`
- **`folds`**: Used by `tune_sparse_keepX()` and `tune_sparse_o2pls()`

### Parallel Processing
- **`ncores`**: Used by `o2cv()` and `o2cv_sparse()`
- **None**: `tune_sparse_keepX()` and `tune_sparse_o2pls()` don't support parallel processing
- **NOT `nr_cores`**: This parameter doesn't exist in any function

---

## Common Mistakes to Avoid

### ❌ Wrong
```r
# Don't use nr_cores (doesn't exist)
tune_sparse_o2pls(..., nr_cores = 4)

# Don't use nr_folds with tuning functions
tune_sparse_keepX(..., nr_folds = 5)

# Don't use folds with stability_selection
stability_selection(..., folds = 10)
```

### ✅ Correct
```r
# Use ncores for CV functions
o2cv(..., ncores = 4)
o2cv_sparse(..., ncores = 4)

# Use folds (no 'nr_') for tuning functions
tune_sparse_keepX(..., folds = 5)
tune_sparse_o2pls(..., folds = 5)

# Use parallel/ncores for stability_selection
stability_selection(..., parallel = TRUE, ncores = 4)
```

---

## Plotting Differences

### tune_sparse_keepX() - Requires Custom Functions
```r
# Returns simple list - CANNOT use plot() directly
tune_result <- tune_sparse_keepX(X, group, nc = 2, test.keepX = c(10, 20, 30))

# Must use custom functions
source("plot_tune_results.R")
plot_tune_keepX_results(tune_result, c(10, 20, 30))
```

### tune_sparse_o2pls() - Works with plot()
```r
# Returns TuneResult S4 object - CAN use plot() directly
tune_o2pls <- tune_sparse_o2pls(X, Y, nc_range = 1:2, keepX_range = c(10, 20))

# plot() works!
plot(tune_o2pls, type = "heatmap")  # ✅ Works
plot(tune_o2pls, type = "line")     # ✅ Works
```

---

## Files Updated

- ✅ `comprehensive_tutorial.R` - Fixed incorrect parameters
- ✅ `PLOTTING_TUNING_RESULTS.md` - Added parameter reference table
- ✅ `PARAMETER_CORRECTIONS.md` - This file

---

## Quick Check

To verify correct parameter names for any function:

```r
library(o2plsda)

# Check function signature
args(o2cv)
args(o2cv_sparse)
args(tune_sparse_keepX)
args(tune_sparse_o2pls)
```

---

**Last Updated:** December 1, 2025
**Status:** ✅ All parameters corrected
