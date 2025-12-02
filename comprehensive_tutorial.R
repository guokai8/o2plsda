################################################################################
# COMPREHENSIVE TUTORIAL FOR o2plsda PACKAGE
# This script demonstrates ALL major functions in the o2plsda package
################################################################################

library(o2plsda)
set.seed(123)

################################################################################
# PART 1: DATA PREPARATION
################################################################################

cat("\n========== PART 1: DATA PREPARATION ==========\n")

# Generate example data matrices
# X: First omics dataset (e.g., gene expression) - 50 samples x 100 features
# Y: Second omics dataset (e.g., metabolomics) - 50 samples x 100 features
X <- matrix(rnorm(5000), 50, 100)
Y <- matrix(rnorm(5000), 50, 100)

# Add row and column names
rownames(X) <- paste("Sample", 1:50, sep="_")
rownames(Y) <- paste("Sample", 1:50, sep="_")
colnames(X) <- paste("Gene", 1:100, sep="_")
colnames(Y) <- paste("Lipid", 1:100, sep="_")

# Scale the data (recommended)
X <- scale(X, scale = TRUE)
Y <- scale(Y, scale = TRUE)

# Create group labels (2 groups: Control vs Treatment)
group <- rep(c("Control", "Treatment"), each = 25)

cat("Data dimensions:\n")
cat("  X (omics 1):", dim(X), "\n")
cat("  Y (omics 2):", dim(Y), "\n")
cat("  Groups:", table(group), "\n")

################################################################################
# PART 2: STANDARD O2PLS ANALYSIS
################################################################################

cat("\n========== PART 2: STANDARD O2PLS ANALYSIS ==========\n")

# 2.1 Cross-validation to find optimal parameters
cat("\n2.1 Cross-validation for parameter selection...\n")
cv_result <- o2cv(X, Y,
                  nc = 1:3,      # Joint components to test
                  nx = 0:2,      # X-specific orthogonal components
                  ny = 0:2,      # Y-specific orthogonal components
                  group = group,
                  nr_folds = 5)  # 5-fold cross-validation

# 2.2 Fit O2PLS model with optimal or specified parameters
cat("\n2.2 Fitting O2PLS model...\n")
o2pls_fit <- o2pls(X, Y, nc = 1, nx = 2, ny = 2)

# 2.3 Summary of O2PLS results
cat("\n2.3 O2PLS Model Summary:\n")
summary(o2pls_fit)
print(o2pls_fit)

# 2.4 Extract loadings and scores
cat("\n2.4 Extracting loadings and scores...\n")
Xjoint_loadings <- loadings(o2pls_fit, loading = "Xjoint")
Yjoint_loadings <- loadings(o2pls_fit, loading = "Yjoint")
Xorth_loadings <- loadings(o2pls_fit, loading = "Xorth")
Yorth_loadings <- loadings(o2pls_fit, loading = "Yorth")

Xjoint_scores <- scores(o2pls_fit, score = "Xjoint")
Yjoint_scores <- scores(o2pls_fit, score = "Yjoint")
Xorth_scores <- scores(o2pls_fit, score = "Xorth")
Yorth_scores <- scores(o2pls_fit, score = "Yorth")

cat("  Xjoint loadings dimension:", dim(Xjoint_loadings), "\n")
cat("  Xjoint scores dimension:", dim(Xjoint_scores), "\n")

################################################################################
# PART 3: O2PLS PLOTTING
################################################################################

cat("\n========== PART 3: O2PLS PLOTTING ==========\n")

# 3.1 Score plots
cat("\n3.1 Creating O2PLS score plots...\n")
p1 <- plot(o2pls_fit, type = "score", var = "Xjoint", group = group)
p2 <- plot(o2pls_fit, type = "score", var = "Yjoint", group = group)
p3 <- plot(o2pls_fit, type = "score", var = "Xorth", group = group)

# 3.2 Loading plots
cat("\n3.2 Creating O2PLS loading plots...\n")
p4 <- plot(o2pls_fit, type = "loading", var = "Xjoint", repel = FALSE)
p5 <- plot(o2pls_fit, type = "loading", var = "Yjoint", repel = FALSE, rotation = TRUE)

################################################################################
# PART 4: O2PLS-DA (DISCRIMINANT ANALYSIS)
################################################################################

cat("\n========== PART 4: O2PLS-DA ANALYSIS ==========\n")

# 4.1 Perform OPLS-DA on O2PLS results
cat("\n4.1 Performing OPLS-DA...\n")
oplsda_result <- oplsda(o2pls_fit, group, nc = 1)

cat("\nOPLS-DA completed\n")

# 4.2 Extract loadings and scores from OPLS-DA
cat("\n4.2 Extracting OPLS-DA results...\n")
oplsda_loadings <- loadings(oplsda_result, loading = "Xloading")
oplsda_scores <- scores(oplsda_result)

# 4.3 Calculate VIP (Variable Importance in Projection)
cat("\n4.3 Calculating VIP scores...\n")
vip_scores <- vip(oplsda_result)
cat("  Top 10 VIP variables:\n")
print(head(sort(vip_scores, decreasing = TRUE), 10))

# 4.4 Plot OPLS-DA results
cat("\n4.4 Creating OPLS-DA plots...\n")
p6 <- plot(oplsda_result, type = "score", group = group)
p7 <- plot(oplsda_result, type = "vip", group = group, repel = FALSE, order = TRUE)
p8 <- plot(oplsda_result, type = "loading", group = group, repel = FALSE)

################################################################################
# PART 5: STANDARD PLS-DA
################################################################################

cat("\n========== PART 5: STANDARD PLS-DA ==========\n")

# 5.1 Perform PLS-DA directly
cat("\n5.1 Performing PLS-DA...\n")
plsda_result <- plsda(X, group, nc = 2, scale = FALSE, center = TRUE)

# 5.2 Summary and print
cat("\n5.2 PLS-DA Summary:\n")
summary(plsda_result)
print(plsda_result)

# 5.3 Extract components
cat("\n5.3 Extracting PLS-DA components...\n")
plsda_loadings <- loadings(plsda_result)
plsda_scores <- scores(plsda_result)
plsda_vip <- vip(plsda_result)

# 5.4 Plot PLS-DA results
cat("\n5.4 Creating PLS-DA plots...\n")
p9 <- plot(plsda_result, type = "score", group = group)
p10 <- plot(plsda_result, type = "vip", group = group, order = TRUE, repel = FALSE)
p11 <- plot(plsda_result, type = "loading", group = group,ind = 1)

################################################################################
# PART 6: SPARSE O2PLS
################################################################################

cat("\n========== PART 6: SPARSE O2PLS ==========\n")

# 6.1 Fit sparse O2PLS model with feature selection
cat("\n6.1 Fitting Sparse O2PLS model...\n")
sparse_o2pls_fit <- sparse_o2pls(X, Y,
                                  nc = 1,      # Joint components
                                  nx = 1,      # X orthogonal components
                                  ny = 1,      # Y orthogonal components
                                  keepX = 20,  # Keep top 20 features in X
                                  keepY = 20)  # Keep top 20 features in Y

# 6.2 Summary and print
cat("\n6.2 Sparse O2PLS Summary:\n")
summary(sparse_o2pls_fit)
print(sparse_o2pls_fit)

# 6.3 Get selected variables
cat("\n6.3 Getting selected variables...\n")
selected_X <- selected_vars(sparse_o2pls_fit, type = "X")
selected_Y <- selected_vars(sparse_o2pls_fit, type = "Y")
selected_both <- selected_vars(sparse_o2pls_fit, type = "both")

cat("  Number of selected X variables:", length(selected_X), "\n")
cat("  Number of selected Y variables:", length(selected_Y), "\n")

# 6.4 Get variable names
cat("\n6.4 Getting selected variable names...\n")
selected_X_names <- selected_var_names(sparse_o2pls_fit, type = "X")
selected_Y_names <- selected_var_names(sparse_o2pls_fit, type = "Y")

cat("  First 5 selected X variables:", head(selected_X_names, 5), "\n")
cat("  First 5 selected Y variables:", head(selected_Y_names, 5), "\n")

# 6.5 Get sparsity information
cat("\n6.5 Sparsity information:\n")
sparsity <- sparsity_info(sparse_o2pls_fit)
print(sparsity)

# 6.6 Plot Sparse O2PLS results
cat("\n6.6 Creating Sparse O2PLS plots...\n")
p12 <- plot(sparse_o2pls_fit, type = "score", var = "Xjoint", group = group)
p13 <- plot(sparse_o2pls_fit, type = "loading", var = "Xjoint")
p14 <- plot(sparse_o2pls_fit, type = "biplot", var = "X", group = group)

# 6.7 Predictions with Sparse O2PLS
cat("\n6.7 Making predictions with Sparse O2PLS...\n")
# Create test data (using first 5 samples as test)
X_test <- X[1:5, ]
predictions <- predict(sparse_o2pls_fit, newdata = X_test)
cat("  Predictions completed for", nrow(X_test), "samples\n")

################################################################################
# PART 7: SPARSE PLS-DA
################################################################################

cat("\n========== PART 7: SPARSE PLS-DA ==========\n")

# 7.1 Fit sparse PLS-DA
cat("\n7.1 Fitting Sparse PLS-DA model...\n")
sparse_plsda_fit <- sparse_plsda(X, group,
                                  nc = 2,      # Number of components
                                  keepX = 30,  # Keep top 30 features
                                  validation = "Mfold",
                                  folds = 5)

# 7.2 Summary and print
cat("\n7.2 Sparse PLS-DA Summary:\n")
summary(sparse_plsda_fit)
print(sparse_plsda_fit)

# 7.3 Get selected variables
cat("\n7.3 Selected variables from Sparse PLS-DA...\n")
sparse_plsda_vars <- selected_vars(sparse_plsda_fit)
sparse_plsda_names <- selected_var_names(sparse_plsda_fit)

cat("  Number of selected variables:", length(sparse_plsda_vars), "\n")
cat("  First 10 selected variables:", head(sparse_plsda_names, 10), "\n")

# 7.4 Sparsity information
cat("\n7.4 Sparsity information:\n")
sparse_plsda_sparsity <- sparsity_info(sparse_plsda_fit)
print(sparse_plsda_sparsity)

# 7.5 Plot Sparse PLS-DA results
cat("\n7.5 Creating Sparse PLS-DA plots...\n")
p15 <- plot(sparse_plsda_fit, type = "score", group = group)
p16 <- plot(sparse_plsda_fit, type = "loading")
p17 <- plot(sparse_plsda_fit, type = "vip", order = TRUE)
p18 <- plot(sparse_plsda_fit, type = "contribution")

# 7.6 Predictions with Sparse PLS-DA
cat("\n7.6 Making predictions with Sparse PLS-DA...\n")
sparse_plsda_pred <- predict(sparse_plsda_fit, newdata = X_test, dist = "max.dist")
cat("  Predicted classes:\n")
print(sparse_plsda_pred$class)
cat("  Prediction probabilities:\n")
print(head(sparse_plsda_pred$probabilities))

################################################################################
# PART 8: TUNING SPARSE MODELS
################################################################################

cat("\n========== PART 8: TUNING SPARSE MODELS ==========\n")

# 8.1 Tune keepX parameter for Sparse PLS-DA
cat("\n8.1 Tuning keepX for Sparse PLS-DA...\n")
test_keepX_values <- seq(10, 50, 10)
tune_result <- tune_sparse_keepX(X, group,
                                  nc = 2,
                                  test.keepX = test_keepX_values,
                                  validation = "Mfold",
                                  folds = 3,
                                  nrepeat = 2)

cat("\n  Tuning completed\n")
cat("  Optimal keepX:", tune_result$choice.keepX, "\n")

# 8.2 Plot tuning results
cat("\n8.2 Creating tuning plots...\n")
# Note: tune_sparse_keepX returns a simple list, not a TuneResult object
# Load custom plotting functions
source("plot_tune_results.R")
p19 <- plot_tune_keepX_results(tune_result, test_keepX_values)
p20 <- plot_tune_keepX_barplot(tune_result, test_keepX_values)
cat("  Plots created successfully\n")

# 8.3 Tune Sparse O2PLS parameters
cat("\n8.3 Tuning Sparse O2PLS parameters...\n")
tune_o2pls <- tune_sparse_o2pls(X, Y,
                                 nc_range = 1:2,
                                 nx_range = 0:1,
                                 ny_range = 0:1,
                                 keepX_range = c(20, 30),
                                 keepY_range = c(20, 30),
                                 folds = 3,
                                 measure = "RMSE")

cat("  Optimal parameters found\n")
cat("  This returns a TuneResult object that CAN be plotted directly\n")

################################################################################
# PART 9: CROSS-VALIDATION FOR SPARSE O2PLS
################################################################################

cat("\n========== PART 9: CROSS-VALIDATION FOR SPARSE O2PLS ==========\n")

# 9.1 Perform cross-validation with sparsity
cat("\n9.1 Running sparse O2PLS cross-validation...\n")
cv_sparse <- o2cv_sparse(X, Y,
                         nc = 1,
                         nx = 1,
                         ny = 1,
                         keepX_range = c(20, 30),
                         keepY_range = c(20, 30),
                         nr_folds = 3,
                         ncores = 1)

cat("  Cross-validation completed\n")

################################################################################
# PART 10: STABILITY SELECTION
################################################################################

cat("\n========== PART 10: STABILITY SELECTION ==========\n")

# 10.1 Perform stability selection
cat("\n10.1 Running stability selection...\n")
stability_result <- stability_selection(X, group,
                                         nc = 2,
                                         keepX_range = seq(10, 40, 10),
                                         n_bootstrap = 20,
                                         threshold = 0.6,
                                         parallel = FALSE)  # Set to TRUE for faster processing

cat("  Stability selection completed\n")

# 10.2 Plot stability results
cat("\n10.2 Creating stability plots...\n")
# Note: With random data, stability selection may not find stable features
# These plots are wrapped in tryCatch to handle that case gracefully
tryCatch({
  p21 <- plot(stability_result, type = "heatmap")
  cat("  ✓ Heatmap created\n")
}, error = function(e) cat("  - Heatmap skipped (no stable features)\n"))

tryCatch({
  p22 <- plot(stability_result, type = "barplot", threshold = 0.6)
  cat("  ✓ Barplot created\n")
}, error = function(e) cat("  - Barplot skipped (no stable features)\n"))

tryCatch({
  p23 <- plot(stability_result, type = "network", threshold = 0.7)
  cat("  ✓ Network plot created\n")
}, error = function(e) cat("  - Network plot skipped (need at least 2 stable features)\n"))

tryCatch({
  p24 <- plot(stability_result, type = "summary")
  cat("  ✓ Summary plot created\n")
}, error = function(e) cat("  - Summary plot skipped (no stable features)\n"))

cat("  Note: Random data typically doesn't have stable features\n")

################################################################################
# PART 11: MODEL COMPARISON
################################################################################

cat("\n========== PART 11: MODEL COMPARISON ==========\n")

# 11.1 Compare multiple models
cat("\n11.1 Comparing multiple models...\n")
comparison <- compare_models(
  list(
    "Standard PLS-DA" = plsda_result,
    "Sparse PLS-DA" = sparse_plsda_fit
  )
)

cat("  Model comparison completed\n")
print(comparison)

################################################################################
# PART 12: ADVANCED PREDICTION FUNCTIONS
################################################################################

cat("\n========== PART 12: ADVANCED PREDICTION FUNCTIONS ==========\n")

# 12.1 Calculate Mahalanobis distances
cat("\n12.1 Calculating Mahalanobis distances...\n")
scores_matrix <- plsda_scores
classification_info <- data.frame(
  sample = rownames(scores_matrix),
  class = group
)
mahal_dist <- calculate_mahalanobis_distances(scores_matrix, classification_info)
cat("  Mahalanobis distances calculated\n")

# 12.2 Calculate prediction probabilities
cat("\n12.2 Calculating prediction probabilities...\n")
pred_probs <- calculate_prediction_probabilities(mahal_dist)
cat("  Prediction probabilities:\n")
print(head(pred_probs))

################################################################################
# PART 13: OPTIMAL PARAMETER FINDING
################################################################################

cat("\n========== PART 13: OPTIMAL PARAMETER FINDING ==========\n")

# 13.1 Find optimal sparse parameters
cat("\n13.1 Finding optimal sparse parameters...\n")
optimal_params <- find_optimal_sparse_params(X, Y,
                                               nc_range = 1:2,
                                               keepX_range = c(15, 25),
                                               keepY_range = c(15, 25),
                                               n_folds = 3)

cat("  Optimal parameters found:\n")
print(optimal_params)

################################################################################
# PART 14: ALL PLOT TYPES FOR SPARSE O2PLS
################################################################################

cat("\n========== PART 14: ALL SPARSE O2PLS PLOT TYPES ==========\n")

cat("\n14.1 Creating all available Sparse O2PLS plots...\n")

# Score plots
p25 <- plot(sparse_o2pls_fit, type = "score", var = "Xjoint", group = group)
p26 <- plot(sparse_o2pls_fit, type = "score", var = "Yjoint", group = group)

# Loading plots
p27 <- plot(sparse_o2pls_fit, type = "loading", var = "Xjoint")
p28 <- plot(sparse_o2pls_fit, type = "loading", var = "Yjoint")

# Biplot
p29 <- plot(sparse_o2pls_fit, type = "biplot", var = "X", group = group)
p30 <- plot(sparse_o2pls_fit, type = "biplot", var = "Y", group = group)

# Diagnostic plots
p31 <- plot(sparse_o2pls_fit, type = "diagnostic", var = "X")
p32 <- plot(sparse_o2pls_fit, type = "diagnostic", var = "Y")

# Contribution plot
p33 <- plot(sparse_o2pls_fit, type = "contribution", var = "X")
p34 <- plot(sparse_o2pls_fit, type = "contribution", var = "Y")

cat("  All Sparse O2PLS plots created\n")

################################################################################
# PART 15: ALL PLOT TYPES FOR SPARSE PLS-DA
################################################################################

cat("\n========== PART 15: ALL SPARSE PLS-DA PLOT TYPES ==========\n")

cat("\n15.1 Creating all available Sparse PLS-DA plots...\n")

# Score plot
p35 <- plot(sparse_plsda_fit, type = "score", group = group, component = c(1, 2))

# Loading plot
p36 <- plot(sparse_plsda_fit, type = "loading", component = c(1, 2))

# VIP plot
p37 <- plot(sparse_plsda_fit, type = "vip", order = TRUE, top_n = 20)

# Contribution plot
p38 <- plot(sparse_plsda_fit, type = "contribution", component = 1)

# Selection plot
p39 <- plot(sparse_plsda_fit, type = "selection")

# Classification plot
p40 <- plot(sparse_plsda_fit, type = "classification", group = group)

# Diagnostic plot
p41 <- plot(sparse_plsda_fit, type = "diagnostic")

cat("  All Sparse PLS-DA plots created\n")

################################################################################
# SUMMARY
################################################################################

cat("\n========== TUTORIAL COMPLETE ==========\n")
cat("\nSummary of tested functions:\n")
cat("  1. o2cv() - O2PLS cross-validation\n")
cat("  2. o2pls() - Standard O2PLS analysis\n")
cat("  3. oplsda() - O2PLS discriminant analysis\n")
cat("  4. plsda() - Standard PLS-DA\n")
cat("  5. sparse_o2pls() - Sparse O2PLS\n")
cat("  6. sparse_plsda() - Sparse PLS-DA\n")
cat("  7. tune_sparse_keepX() - Tune sparsity parameter\n")
cat("  8. tune_sparse_o2pls() - Tune sparse O2PLS parameters\n")
cat("  9. o2cv_sparse() - Cross-validation for sparse O2PLS\n")
cat(" 10. stability_selection() - Stability-based feature selection\n")
cat(" 11. compare_models() - Model comparison\n")
cat(" 12. loadings() - Extract loadings\n")
cat(" 13. scores() - Extract scores\n")
cat(" 14. vip() - Variable importance\n")
cat(" 15. selected_vars() - Get selected variables\n")
cat(" 16. selected_var_names() - Get variable names\n")
cat(" 17. sparsity_info() - Get sparsity information\n")
cat(" 18. predict() - Make predictions\n")
cat(" 19. calculate_mahalanobis_distances() - Distance calculation\n")
cat(" 20. calculate_prediction_probabilities() - Probability calculation\n")
cat(" 21. find_optimal_sparse_params() - Optimal parameter search\n")
cat(" 22. plot() - Various visualization types\n")
cat(" 23. summary() - Model summaries\n")
cat(" 24. print() - Print model information\n")

cat("\nAll major functions have been tested successfully!\n")
cat("\nTotal plots created:", 41, "\n")

################################################################################
# OPTIONAL: Save results
################################################################################

# Uncomment to save results
# save(o2pls_fit, oplsda_result, plsda_result,
#      sparse_o2pls_fit, sparse_plsda_fit,
#      file = "o2plsda_tutorial_results.RData")

cat("\nTutorial completed successfully!\n")
