# Testing Guide for o2plsda Package

## Files Created

Three testing files have been created for you:

### 1. **quick_test.R** - Fast Verification (2-3 minutes)
Quick sanity check to verify all major functions work.

**Run with:**
```r
source("quick_test.R")
```

**What it tests:**
- ✓ o2pls() - O2PLS analysis
- ✓ o2cv() - Cross-validation
- ✓ oplsda() - OPLS-DA
- ✓ plsda() - PLS-DA
- ✓ sparse_o2pls() - Sparse O2PLS
- ✓ sparse_plsda() - Sparse PLS-DA
- ✓ predict() - Predictions
- ✓ tune_sparse_keepX() - Tuning
- ✓ stability_selection() - Stability
- ✓ All plot types

**Expected output:** Progress messages with checkmarks (✓) for each successful test.

---

### 2. **comprehensive_tutorial.R** - Complete Testing (10-15 minutes)
Comprehensive tutorial testing ALL package functions with detailed examples.

**Run with:**
```r
source("comprehensive_tutorial.R")
```

**What it covers (15 parts):**
1. Data Preparation
2. Standard O2PLS Analysis
3. O2PLS Plotting
4. O2PLS-DA (Discriminant Analysis)
5. Standard PLS-DA
6. Sparse O2PLS
7. Sparse PLS-DA
8. Tuning Sparse Models
9. Cross-Validation for Sparse O2PLS
10. Stability Selection
11. Model Comparison
12. Advanced Prediction Functions
13. Optimal Parameter Finding
14. All Sparse O2PLS Plot Types
15. All Sparse PLS-DA Plot Types

**Creates:** 41+ plot objects (p1, p2, ..., p41)

**Tests:** 24 main functions with all their parameters

---

### 3. **TUTORIAL_GUIDE.md** - Documentation
Detailed documentation explaining every aspect of the tutorial.

**Contains:**
- Overview of all 15 tutorial parts
- Function descriptions
- Plot type references
- Customization examples
- Troubleshooting tips
- Performance optimization advice
- How to adapt for your own data

---

## Quick Start

### Option A: Just verify everything works
```r
# Run quick test (2-3 minutes)
source("quick_test.R")

# If all tests pass with ✓, your package is working!
```

### Option B: Learn and test everything
```r
# Run comprehensive tutorial (10-15 minutes)
source("comprehensive_tutorial.R")

# View created plots
print(p1)   # First plot
print(p15)  # Sparse PLS-DA plot
# ... up to p41
```

### Option C: Run specific sections
```r
# Open comprehensive_tutorial.R in RStudio
# Run sections individually with Ctrl+Enter or Cmd+Enter
# Start with PART 1 (Data Preparation)
# Then run any other parts you're interested in
```

---

## Function Coverage

### Core Analysis Functions (5)
| Function | Tested in Quick | Tested in Comprehensive | Description |
|----------|-----------------|-------------------------|-------------|
| `o2pls()` | ✓ | ✓ | Standard O2PLS |
| `oplsda()` | ✓ | ✓ | O2PLS-DA |
| `plsda()` | ✓ | ✓ | PLS-DA |
| `sparse_o2pls()` | ✓ | ✓ | Sparse O2PLS |
| `sparse_plsda()` | ✓ | ✓ | Sparse PLS-DA |

### Cross-Validation & Tuning (5)
| Function | Tested in Quick | Tested in Comprehensive | Description |
|----------|-----------------|-------------------------|-------------|
| `o2cv()` | ✓ | ✓ | O2PLS CV |
| `o2cv_sparse()` | - | ✓ | Sparse O2PLS CV |
| `tune_sparse_keepX()` | ✓ | ✓ | Tune sparsity |
| `tune_sparse_o2pls()` | - | ✓ | Tune all params |
| `find_optimal_sparse_params()` | - | ✓ | Auto optimization |

### Feature Selection (4)
| Function | Tested in Quick | Tested in Comprehensive | Description |
|----------|-----------------|-------------------------|-------------|
| `stability_selection()` | ✓ | ✓ | Stability selection |
| `selected_vars()` | ✓ | ✓ | Get indices |
| `selected_var_names()` | - | ✓ | Get names |
| `sparsity_info()` | ✓ | ✓ | Sparsity stats |

### Extraction Functions (3)
| Function | Tested in Quick | Tested in Comprehensive | Description |
|----------|-----------------|-------------------------|-------------|
| `loadings()` | ✓ | ✓ | Extract loadings |
| `scores()` | ✓ | ✓ | Extract scores |
| `vip()` | ✓ | ✓ | Variable importance |

### Prediction (3)
| Function | Tested in Quick | Tested in Comprehensive | Description |
|----------|-----------------|-------------------------|-------------|
| `predict()` | ✓ | ✓ | Make predictions |
| `calculate_mahalanobis_distances()` | - | ✓ | Distance calc |
| `calculate_prediction_probabilities()` | - | ✓ | Probability calc |

### Utilities (4)
| Function | Tested in Quick | Tested in Comprehensive | Description |
|----------|-----------------|-------------------------|-------------|
| `compare_models()` | - | ✓ | Model comparison |
| `summary()` | ✓ | ✓ | Model summaries |
| `print()` | Implicit | ✓ | Print info |
| `plot()` | ✓ | ✓ | Visualization |

**Total: 24 functions comprehensively tested**

---

## Plot Types Tested

### O2PLS Plots (5 types)
- Score plots (Xjoint, Yjoint, Xorth, Yorth)
- Loading plots (with rotation, repel options)

### O2PLS-DA Plots (3 types)
- Score plot
- VIP plot
- Loading plot

### PLS-DA Plots (3 types)
- Score plot
- VIP plot
- Loading plot

### Sparse O2PLS Plots (10+ types)
- Score plots (Xjoint, Yjoint)
- Loading plots (Xjoint, Yjoint)
- Biplots (X, Y)
- Diagnostic plots (X, Y)
- Contribution plots (X, Y)

### Sparse PLS-DA Plots (7 types)
- Score plot
- Loading plot
- VIP plot
- Contribution plot
- Selection plot
- Classification plot
- Diagnostic plot

### Meta-Analysis Plots (4 types)
- Tuning heatmap
- Tuning line plot
- Stability heatmap
- Stability barplot
- Stability network
- Stability summary

**Total: 30+ different plot types tested**

---

## Validation Checklist

Use this checklist to verify your testing:

### Basic Functionality
- [ ] quick_test.R runs without errors
- [ ] All 10 test sections pass with ✓
- [ ] Plots can be displayed with `print(p1)`, etc.

### Core Methods
- [ ] O2PLS fits successfully
- [ ] Cross-validation completes
- [ ] OPLS-DA produces results
- [ ] PLS-DA works
- [ ] VIP scores are calculated

### Sparse Methods
- [ ] Sparse O2PLS selects features
- [ ] Sparse PLS-DA selects features
- [ ] Selected variables can be extracted
- [ ] Sparsity information is available

### Predictions
- [ ] Sparse O2PLS predictions work
- [ ] Sparse PLS-DA predictions work
- [ ] Predicted classes are reasonable
- [ ] Probabilities are calculated

### Tuning & Selection
- [ ] keepX tuning completes
- [ ] Stability selection runs
- [ ] Results can be plotted

### Visualizations
- [ ] Score plots render correctly
- [ ] Loading plots show features
- [ ] VIP plots display importance
- [ ] Heatmaps show patterns

---

## Expected Runtime

**Hardware specs for estimates:** Standard laptop with 4 cores

| Script | Expected Time | Memory Usage |
|--------|---------------|--------------|
| quick_test.R | 2-3 minutes | ~200 MB |
| comprehensive_tutorial.R | 10-15 minutes | ~500 MB |
| Individual sections | 30s - 2 min each | ~200 MB |

**Note:** Times may vary based on:
- Number of CPU cores available
- Whether parallel processing is used
- System memory
- R version and package versions

---

## Troubleshooting

### If quick_test.R fails:

1. **Check package installation:**
   ```r
   library(o2plsda)
   packageVersion("o2plsda")
   ```

2. **Check dependencies:**
   ```r
   # Required packages
   library(ggplot2)
   library(dplyr)
   library(Rcpp)
   ```

3. **Update packages:**
   ```r
   update.packages(ask = FALSE)
   ```

### If plots don't display:

```r
# Explicitly print plot
print(p1)

# Or save to file
ggsave("test_plot.pdf", p1)

# Check plot object
class(p1)  # Should be "gg" "ggplot"
```

### If memory issues occur:

```r
# Run with smaller data
X <- scale(matrix(rnorm(1000), 50, 20))  # Smaller features
Y <- scale(matrix(rnorm(1000), 50, 20))

# Use fewer cross-validation folds
cv <- o2cv(X, Y, 1:2, 0:1, 0:1, nr_folds = 3)  # Instead of 10
```

### If slow performance:

```r
# Enable parallel processing
library(parallel)
ncores <- detectCores() - 1

cv <- o2cv(X, Y, 1:3, 0:2, 0:2, nr_cores = ncores)
```

---

## Next Steps After Testing

1. **Verify all tests pass** ✓
2. **Examine some plots** to understand output
3. **Read TUTORIAL_GUIDE.md** for detailed explanations
4. **Adapt code for your own data**
5. **Customize plots** for publications
6. **Export results** as needed

---

## Saving Results

### Save all analysis results:
```r
save(o2pls_fit, oplsda_result, plsda_result,
     sparse_o2pls_fit, sparse_plsda_fit,
     file = "analysis_results.RData")
```

### Save selected variables:
```r
selected_features <- data.frame(
  variable = selected_var_names(sparse_plsda_fit),
  index = selected_vars(sparse_plsda_fit)
)
write.csv(selected_features, "selected_features.csv", row.names = FALSE)
```

### Save all plots:
```r
# To PDF
pdf("all_test_plots.pdf", width = 10, height = 8)
for(i in 1:41) {
  plot_obj <- get(paste0("p", i))
  print(plot_obj)
}
dev.off()

# Individual plots
ggsave("o2pls_scores.pdf", p1, width = 8, height = 6)
ggsave("sparse_plsda_scores.pdf", p15, width = 8, height = 6)
```

---

## Support

If you encounter issues:

1. Check this troubleshooting guide
2. Review TUTORIAL_GUIDE.md
3. Check function help: `?o2pls`, `?sparse_plsda`, etc.
4. Visit: https://github.com/guokai8/o2plsda/issues
5. Email: guokai8@gmail.com

---

## Package Information

- **Package:** o2plsda
- **Version:** 0.0.28
- **Maintainer:** Kai Guo
- **License:** GPL-3
- **URL:** https://github.com/guokai8/o2plsda

---

**Testing prepared:** December 2025

**Status:** Ready for comprehensive testing ✓
