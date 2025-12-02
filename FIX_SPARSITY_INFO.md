# Fix for sparsity_info() Error with sparse_plsda Objects

## The Problem

When calling `sparsity_info()` on a `sparse_plsda` object:

```r
sparse_plsda_fit <- sparse_plsda(X, group, nc = 2, keepX = 30)
sparse_plsda_sparsity <- sparsity_info(sparse_plsda_fit)
```

**Error:**
```
Error in x$results : Unknown slot/element: results
```

## Root Cause

Two issues were causing this error:

### Issue 1: Conflicting Internal Function
There was an internal `sparsity_info()` function in `plot.R` (line 1294) that conflicted with the exported S3 generic in `AllGenerics.R`. This internal function was trying to access `x$results$Xloading` which doesn't exist for sparse_plsda objects.

### Issue 2: S3/S4 Method Dispatch
`sparsity_info` was defined as an S3 generic, but `sparse_plsda` is an S4 class. S3 method dispatch doesn't work automatically with S4 classes.

## The Fixes

### Fix 1: Renamed Internal Function
Renamed the conflicting internal function in `plot.R`:
```r
# Before (line 1294)
sparsity_info <- function(x) { ... }

# After
get_sparsity_info_internal <- function(x) { ... }
```

### Fix 2: Updated sparsity_info.sparse_plsda()
Fixed the S3 method to use direct slot access (`@`) for S4 objects:
```r
sparsity_info.sparse_plsda <- function(x, ...) {
    if(inherits(x, "sparse_plsda")) {
        # S4 object - use @ to access slots
        total_vars <- ncol(x@X)
        selected_vars <- length(x@selected_vars)
        sparsity <- 1 - selected_vars / total_vars
        keepX_val <- x@keepX
    } else {
        # List object (legacy)
        total_vars <- if(!is.null(x$total_variables)) x$total_variables else nrow(x$loadings)
        selected_vars <- length(x$selected_vars)
        sparsity <- 1 - selected_vars / total_vars
        keepX_val <- x$keepX
    }

    return(list(
        selected = selected_vars,
        total = total_vars,
        sparsity = sparsity,
        percent_kept = round((1 - sparsity) * 100, 1),
        keepX = keepX_val
    ))
}
```

### Fix 3: Added S4 Method
Added an S4 method to properly handle S4 class dispatch:
```r
setGeneric("sparsity_info")

setMethod("sparsity_info", "sparse_plsda", function(x, ...) {
    sparsity_info.sparse_plsda(x, ...)
})
```

## Now It Works!

```r
library(o2plsda)
set.seed(123)

# Create data
X <- scale(matrix(rnorm(2500), 50, 50))
group <- rep(c("A", "B"), each = 25)

# Fit model
sparse_plsda_fit <- sparse_plsda(X, group, nc = 1, keepX = 30)

# Get sparsity info - NOW WORKS!
sparsity <- sparsity_info(sparse_plsda_fit)
print(sparsity)
```

**Output:**
```
$selected
[1] 30

$total
[1] 50

$sparsity
[1] 0.4

$percent_kept
[1] 60

$keepX
[1] 30
```

## What You Get

The `sparsity_info()` function returns a list with:

- **selected**: Number of variables selected
- **total**: Total number of variables in the original dataset
- **sparsity**: Proportion of variables removed (0 = no sparsity, 1 = all removed)
- **percent_kept**: Percentage of variables kept (100 - sparsity × 100)
- **keepX**: The keepX parameter used in the model

## Related Functions

These functions also work correctly with sparse_plsda objects:

```r
# Get selected variable indices
vars <- selected_vars(sparse_plsda_fit)

# Get selected variable names
var_names <- selected_var_names(sparse_plsda_fit)

# Summary
summary(sparse_plsda_fit)
```

## Files Modified

1. **R/AllGenerics.R**
   - Fixed `sparsity_info.sparse_plsda()` to use `@` for S4 slots
   - Added `setGeneric("sparsity_info")`
   - Added `setMethod("sparsity_info", "sparse_plsda", ...)`

2. **R/plot.R**
   - Renamed `sparsity_info()` to `get_sparsity_info_internal()` (line 1294)

## Status

✅ **FIXED and TESTED**

The package has been rebuilt and all sparsity_info functions work correctly with:
- sparse_plsda (S4 class)
- sparse_o2pls (list)
- SparseO2pls (S4 class)

---

**Fixed:** December 1, 2025
**Package Version:** 0.0.28
