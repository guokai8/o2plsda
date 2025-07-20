#' Sparse Partial Least Squares Discriminant Analysis
#' 
#' Performs sparse PLS-DA for classification with automatic variable selection.
#' 
#' @param X Numeric matrix (samples x variables). The predictor data matrix.
#' @param Y Factor or character vector. Class labels for each sample.
#' @param nc Integer. Number of components to extract.
#' @param keepX Vector of integers. Number of variables to keep per component.
#'   If NULL, automatically determined.
#' @param validation Character. Validation method for parameter tuning: "Mfold" or "loo".
#' @param folds Integer. Number of cross-validation folds (default: 5).
#' @param test.keepX Vector of keepX values to test during tuning.
#' @param scale Logical. Scale variables (default: TRUE).
#' @param center Logical. Center variables (default: TRUE).
#' @param dist Character. Distance metric for classification: "max.dist", 
#'   "centroids.dist", or "mahalanobis.dist".
#' @param tune.keepX Logical. Automatically tune keepX parameters (default: TRUE if keepX is NULL).
#' 
#' @return sparse_plsda object containing sparse loadings, scores, VIP scores, and classification results.
#' 
#' @examples
#' # Example 1: Basic sparse PLS-DA
#' set.seed(123)
#' n <- 120
#' p <- 200  # High-dimensional feature space
#' 
#' # Generate classification data
#' X <- matrix(rnorm(n * p), n, p)
#' 
#' # Create three classes with different patterns
#' n_per_class <- n / 3
#' classes <- rep(c("Class_A", "Class_B", "Class_C"), each = n_per_class)
#' 
#' # Add class-specific signals to subsets of variables
#' X[classes == "Class_A", 1:20] <- X[classes == "Class_A", 1:20] + 2
#' X[classes == "Class_B", 21:40] <- X[classes == "Class_B", 21:40] + 2
#' X[classes == "Class_C", 41:60] <- X[classes == "Class_C", 41:60] + 2
#' 
#' # Add some noise correlation
#' X[, 61:80] <- X[, 1:20] + matrix(rnorm(n * 20, sd = 0.5), n, 20)
#' 
#' colnames(X) <- paste0("Feature_", 1:p)
#' Y <- factor(classes)
#' 
#' # Fit sparse PLS-DA with automatic tuning
#' sparse_plsda_fit <- sparse_plsda(
#'     X = X,
#'     Y = Y,
#'     nc = 2,
#'     validation = "Mfold",
#'     folds = 5,
#'     test.keepX = seq(10, 50, 10),  # Test different sparsity levels
#'     tune.keepX = TRUE
#' )
#' 
#' # View results
#' print(sparse_plsda_fit)
#' summary(sparse_plsda_fit)
#' 
#' # Extract selected features
#' selected_features <- selected_vars(sparse_plsda_fit)
#' selected_names <- selected_var_names(sparse_plsda_fit)
#' 
#' cat("Selected", length(selected_features), "out of", p, "features\n")
#' cat("Selected features:", paste(head(selected_names, 15), collapse = ", "), "\n")
#' 
#' # Plot results
#' plot(sparse_plsda_fit, type = "score", group = Y)
#' plot(sparse_plsda_fit, type = "loading", top = 20)
#' 
#' # Example 2: Manual keepX specification
#' sparse_plsda_manual <- sparse_plsda(
#'     X = X,
#'     Y = Y,
#'     nc = 3,
#'     keepX = c(30, 25, 20),  # Decreasing sparsity per component
#'     tune.keepX = FALSE,
#'     dist = "centroids.dist"
#' )
#' 
#' # Compare classification performance
#' cat("Auto-tuned error rate:", sparse_plsda_fit$classification$overall_error, "\n")
#' cat("Manual keepX error rate:", sparse_plsda_manual$classification$overall_error, "\n")
#' 
#' # Example 3: Cross-validation and prediction
#' # Split data
#' train_idx <- sample(1:n, 0.7 * n)
#' test_idx <- setdiff(1:n, train_idx)
#' 
#' X_train <- X[train_idx, ]
#' Y_train <- Y[train_idx]
#' X_test <- X[test_idx, ]
#' Y_test <- Y[test_idx]
#' 
#' # Train model
#' model_train <- sparse_plsda(X_train, Y_train, nc = 2, keepX = c(25, 20))
#' 
#' # Predict test set
#' predictions <- predict(model_train, X_test, dist = "max.dist")
#' 
#' # Evaluate predictions
#' confusion_matrix <- table(Predicted = predictions$class, Actual = Y_test)
#' accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
#' 
#' print(confusion_matrix)
#' cat("Test set accuracy:", round(accuracy, 3), "\n")
#' 
#' @seealso \code{\link{sparse_o2pls}}, \code{\link{stability_selection}}, \code{\link{tune_sparse_keepX}}
#' @author Kai Guo
#' @export
sparse_plsda <- function(X, Y, nc, keepX = NULL, validation = "Mfold", 
                         folds = 5, test.keepX = seq(5, 50, 5),
                         scale = TRUE, center = TRUE, dist = "max.dist",
                         tune.keepX = is.null(keepX)) {
    
    # Input validation and preprocessing
    X <- as.matrix(X)
    Y <- as.factor(Y)
    n_classes <- nlevels(Y)
    
    # CRITICAL FIX: Ensure minimum components for classification
    if(nc > (n_classes - 1)) {
        nc <- min(nc, n_classes - 1)
        warning("nc reduced to ", nc, " (maximum for ", n_classes, " classes)")
    }
    
    if(tune.keepX && is.null(keepX)) {
        # Automatic keepX tuning
        tune_result <- tune_sparse_keepX(X, Y, nc, test.keepX, validation, folds, dist)
        keepX <- tune_result$choice.keepX
        message("Optimal keepX selected: ", paste(keepX, collapse = ", "))
    } else if(is.null(keepX)) {
        # Use default keepX
        keepX <- auto_determine_keepX(X, nc)
    }
    
    # CRITICAL FIX: Ensure keepX length matches nc
    if(length(keepX) == 1) {
        keepX <- rep(keepX, nc)
    } else if(length(keepX) != nc) {
        keepX <- keepX[1:nc]  # Truncate or pad as needed
        warning("keepX adjusted to match nc = ", nc)
    }
    
    # Convert Y to dummy matrix for PLS
    Y_dummy <- model.matrix(~ . - 1, data = data.frame(class = Y))
    
    # Preprocessing with proper attribute storage
    if(scale || center) {
        X_processed <- scale(X, center = center, scale = scale)
        center_values <- if(center) attr(X_processed, "scaled:center") else NULL
        scale_values <- if(scale) attr(X_processed, "scaled:scale") else NULL
    } else {
        X_processed <- X
        center_values <- NULL
        scale_values <- NULL
    }
    
    # CRITICAL FIX: Enhanced sparse PLS fitting
    sparse_fit <- fit_sparse_pls(X_processed, Y_dummy, nc, keepX)
    
    # Calculate VIP scores
    vip_scores <- calculate_sparse_vip(sparse_fit, X_processed, Y_dummy, nc)
    
    # Calculate classification performance
    class_performance <- evaluate_sparse_classification(
        sparse_fit, X_processed, Y, dist
    )

    result <- new("sparse_plsda",
                  call = match.call(),
                  ncomp = as.numeric(nc),
                  keepX = as.numeric(keepX),
                  scores = sparse_fit$scores,
                  loadings = sparse_fit$loadings,
                  vip = vip_scores,
                  selected_vars = as.numeric(which(rowSums(abs(sparse_fit$loadings)) > 1e-10)),
                  classification = class_performance,
                  preprocessing = list(center = center_values, scale = scale_values),
                  dist = as.character(dist),
                  Y_levels = as.character(levels(Y)),
                  X = X,
                  Y = Y)
    return(result)
}


#' Cross-validation for Sparse PLS-DA
#' 
#' Performs cross-validation to tune keepX parameters for sparse PLS-DA
#' 
#' @param X Numeric matrix (samples x variables)
#' @param Y Factor or character vector of class labels
#' @param nc Number of components
#' @param test.keepX Vector of keepX values to test
#' @param validation Validation method ("Mfold" or "loo")
#' @param folds Number of cross-validation folds
#' @param dist Distance metric for classification
#' @param nrepeat Number of repetitions for robust estimation
#' 
#' @return List with optimal keepX and CV results
#' @export
tune_sparse_keepX <- function(X, Y, nc, test.keepX = seq(5, 50, 5), 
                              validation = "Mfold", folds = 5, 
                              dist = "max.dist", nrepeat = 10) {
    
    X <- as.matrix(X)
    Y <- as.factor(Y)
    
    n_keepX <- length(test.keepX)
    cv_results <- array(NA, dim = c(n_keepX, nc, nrepeat))
    
    cat("Tuning keepX parameters...\n")
    pb <- txtProgressBar(min = 0, max = n_keepX * nc * nrepeat, style = 3)
    counter <- 0
    
    for(rep in 1:nrepeat) {
        for(comp in 1:nc) {
            for(k_idx in 1:n_keepX) {
                counter <- counter + 1
                
                # Create keepX vector for this test
                keepX_test <- rep(test.keepX[k_idx], comp)
                
                # Perform cross-validation
                cv_error <- cv_sparse_plsda(X, Y, comp, keepX_test, 
                                                     validation, folds, dist)
                cv_results[k_idx, comp, rep] <- cv_error
                
                setTxtProgressBar(pb, counter)
            }
        }
    }
    close(pb)
    
    # Find optimal parameters using one standard error rule
    mean_errors <- apply(cv_results, c(1, 2), mean, na.rm = TRUE)
    optimal <- find_optimal_sparse_params(mean_errors, test.keepX)
    
    return(optimal)
}

#' Enhanced cross-validation for sparse PLS-DA
#' @keywords internal
cv_sparse_plsda <- function(X, Y, nc, keepX, validation, folds, dist) {
    n <- nrow(X)
    
    # Create folds
    if(validation == "Mfold") {
        fold_ids <- create_stratified_folds(Y, folds)
    } else if(validation == "loo") {
        fold_ids <- 1:n
        folds <- n
    }
    
    cv_errors <- numeric(folds)
    
    for(fold in 1:folds) {
        test_idx <- which(fold_ids == fold)
        train_idx <- setdiff(1:n, test_idx)
        
        if(length(train_idx) < 3 || length(test_idx) == 0) {
            cv_errors[fold] <- 1.0
            next
        }
        
        X_train <- X[train_idx, , drop = FALSE]
        Y_train <- Y[train_idx]
        X_test <- X[test_idx, , drop = FALSE]
        Y_test <- Y[test_idx]
        
        tryCatch({
            # Fit sparse PLS-DA on training data
            fit <- sparse_plsda(X_train, Y_train, nc, keepX = keepX, 
                                tune.keepX = FALSE)
            
            # Predict on test data
            pred <- predict(fit, X_test, dist = dist)
            
            # Calculate error rate
            cv_errors[fold] <- mean(pred$class != Y_test, na.rm = TRUE)
            
        }, error = function(e) {
            cv_errors[fold] <- 1.0  # Maximum error rate
        })
    }
    
    return(mean(cv_errors, na.rm = TRUE))
}


#' Cross-validation for Sparse O2PLS
#' 
#' Performs cross-validation to evaluate sparse O2PLS models
#' 
#' @param X Numeric matrix (samples x variables)
#' @param Y Numeric matrix (samples x variables)
#' @param nc_range Range of joint components to test
#' @param nx_range Range of X-orthogonal components to test
#' @param ny_range Range of Y-orthogonal components to test
#' @param keepX_range Range of keepX values to test
#' @param keepY_range Range of keepY values to test
#' @param validation Validation method
#' @param folds Number of folds
#' @param measure Performance measure
#' 
#' @return TuneResult object with optimal parameters
#' @export
tune_sparse_o2pls <- function(X, Y, nc_range = 1:3, nx_range = 0:2, ny_range = 0:2,
                              keepX_range = seq(10, 50, 10), keepY_range = seq(10, 50, 10),
                              validation = "Mfold", folds = 5, measure = "RMSE") {
    
    # Create parameter grid including sparsity parameters
    param_grid <- expand.grid(
        nc = nc_range,
        nx = nx_range, 
        ny = ny_range,
        keepX = keepX_range,
        keepY = keepY_range
    )
    
    # Filter valid combinations
    param_grid <- param_grid[param_grid$nc > 0, ]
    param_grid <- param_grid[param_grid$keepX <= ncol(X), ]
    param_grid <- param_grid[param_grid$keepY <= ncol(Y), ]
    
    cat("Testing", nrow(param_grid), "parameter combinations...\n")
    pb <- txtProgressBar(min = 0, max = nrow(param_grid), style = 3)
    
    results <- vector("list", nrow(param_grid))
    
    for(i in 1:nrow(param_grid)) {
        params <- param_grid[i, ]
        
        # Perform cross-validation for this parameter combination
        cv_score <- cv_sparse_o2pls_single(
            X, Y, params$nc, params$nx, params$ny,
            params$keepX, params$keepY, validation, folds, measure
        )
        
        results[[i]] <- list(
            nc = params$nc, nx = params$nx, ny = params$ny,
            keepX = params$keepX, keepY = params$keepY,
            mean_score = cv_score,
            sd_score = NA  # Could be enhanced with multiple repetitions
        )
        
        setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Find optimal parameters
    optimal_params <- find_optimal_parameters(results, measure)
    
    # Create TuneResult object
    tune_result <- new("TuneResult",
                       optimal = optimal_params,
                       all_results = results,
                       method = "sparse_grid_search",
                       measure = measure,
                       validation = validation)
    
    return(tune_result)
}

#' Single cross-validation run for sparse O2PLS
#' @keywords internal
cv_sparse_o2pls_single <- function(X, Y, nc, nx, ny, keepX, keepY, validation, folds, measure) {
    n <- nrow(X)
    
    # Create folds
    if(validation == "Mfold") {
        fold_ids <- sample(rep(1:folds, length.out = n))
    } else if(validation == "loo") {
        fold_ids <- 1:n
        folds <- n
    }
    
    cv_scores <- numeric(folds)
    
    for(fold in 1:folds) {
        test_idx <- which(fold_ids == fold)
        train_idx <- setdiff(1:n, test_idx)
        
        if(length(train_idx) < 3 || length(test_idx) == 0) {
            cv_scores[fold] <- ifelse(measure == "RMSE", Inf, 0)
            next
        }
        
        X_train <- X[train_idx, , drop = FALSE]
        Y_train <- Y[train_idx, , drop = FALSE]
        X_test <- X[test_idx, , drop = FALSE]
        Y_test <- Y[test_idx, , drop = FALSE]
        
        tryCatch({
            # Fit sparse O2PLS
            sparse_fit <- sparse_o2pls(X_train, Y_train, nc, nx, ny,
                                       keepX = rep(keepX, nc + nx),
                                       keepY = rep(keepY, nc + ny))
            
            # Predict Y from X
            Y_pred <- predict(sparse_fit, X_test)
            
            # Calculate error
            if(measure == "RMSE") {
                cv_scores[fold] <- sqrt(mean((Y_test - Y_pred)^2, na.rm = TRUE))
            } else if(measure == "R2") {
                ss_res <- sum((Y_test - Y_pred)^2, na.rm = TRUE)
                ss_tot <- sum((Y_test - mean(Y_train, na.rm = TRUE))^2, na.rm = TRUE)
                cv_scores[fold] <- max(0, 1 - ss_res/ss_tot)
            } else {
                cv_scores[fold] <- mean(abs(Y_test - Y_pred), na.rm = TRUE)
            }
            
        }, error = function(e) {
            cv_scores[fold] <- ifelse(measure %in% c("RMSE", "MAE"), Inf, 0)
        })
    }
    
    return(mean(cv_scores[is.finite(cv_scores)], na.rm = TRUE))
}

#' Enhanced predict method for sparse O2PLS
#' @export
predict.SparseO2pls <- function(object, newdata, ...) {
    if(missing(newdata)) {
        stop("newdata is required for prediction")
    }
    
    newdata <- as.matrix(newdata)
    
    # Apply same preprocessing as training data
    if(!is.null(object@preprocessing$center_X)) {
        newdata <- scale(newdata, center = object@preprocessing$center_X, scale = FALSE)
    }
    if(!is.null(object@preprocessing$scale_X)) {
        newdata <- scale(newdata, center = FALSE, scale = object@preprocessing$scale_X)
    }
    
    # Remove orthogonal components if present
    newdata_proj <- newdata
    
    # Remove X-orthogonal components
    if(ncol(object@results$Xorth_loadings) > 0) {
        for(i in 1:ncol(object@results$Xorth_loadings)) {
            t_orth <- newdata_proj %*% object@results$Xorth_loadings[, i, drop = FALSE]
            newdata_proj <- newdata_proj - t_orth %*% t(object@results$Xorth_loadings[, i, drop = FALSE])
        }
    }
    
    # Calculate joint scores
    X_scores_new <- newdata_proj %*% object@results$Xloading
    
    # Predict Y using regression coefficients
    Y_pred <- X_scores_new %*% object@results$BT %*% t(object@results$Yloading)
    
    # Reverse Y preprocessing if needed
    if(!is.null(object@preprocessing$scale_Y)) {
        Y_pred <- Y_pred * matrix(object@preprocessing$scale_Y, nrow(Y_pred), ncol(Y_pred), byrow = TRUE)
    }
    if(!is.null(object@preprocessing$center_Y)) {
        Y_pred <- Y_pred + matrix(object@preprocessing$center_Y, nrow(Y_pred), ncol(Y_pred), byrow = TRUE)
    }
    
    return(Y_pred)
}

#' Enhanced predict method for sparse PLS-DA
#' @export  
predict.sparse_plsda <- function(object, newdata, dist = "max.dist", ...) {
    if(missing(newdata)) {
        stop("newdata is required for prediction")
    }
    
    newdata <- as.matrix(newdata)
    
    # Apply preprocessing
    if(!is.null(object$preprocessing$center)) {
        newdata <- scale(newdata, center = object$preprocessing$center, scale = FALSE)
    }
    if(!is.null(object$preprocessing$scale)) {
        newdata <- scale(newdata, center = FALSE, scale = object$preprocessing$scale)
    }
    
    # Calculate scores using sparse loadings
    scores_new <- newdata %*% object$loadings
    
    # Classification based on distance metric
    if(dist == "max.dist") {
        distances <- calculate_max_distances(scores_new, object$classification$centroids)
        predicted_classes <- object$Y_levels[apply(distances, 1, which.min)]
    } else if(dist == "centroids.dist") {
        distances <- calculate_centroid_distances(scores_new, object$classification$centroids)
        predicted_classes <- object$Y_levels[apply(distances, 1, which.min)]
    } else if(dist == "mahalanobis.dist") {
        distances <- calculate_mahalanobis_distances(scores_new, object$classification)
        predicted_classes <- object$Y_levels[apply(distances, 1, which.min)]
    } else {
        # Default to centroid distance
        distances <- calculate_centroid_distances(scores_new, object$classification$centroids)
        predicted_classes <- object$Y_levels[apply(distances, 1, which.min)]
    }
    
    # Calculate prediction confidence/probabilities
    prediction_probs <- calculate_prediction_probabilities(distances)
    
    return(list(
        class = predicted_classes,
        scores = scores_new,
        distances = distances,
        probabilities = prediction_probs
    ))
}


#' Enhanced sparse cross-validation that integrates with existing o2cv
#' @export
o2cv_sparse <- function(X, Y, nc, nx, ny, keepX_range = seq(10, 50, 10), 
                        keepY_range = seq(10, 50, 10), group = NULL, 
                        nr_folds = 5, ncores = 1, validation = "Mfold") {
    
    # Input validation
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    
    if(is.null(group)) {
        group <- rep(1, nrow(X))
    }
    
    # Create parameter grid
    sparse_params <- expand.grid(
        keepX = keepX_range,
        keepY = keepY_range
    )
    
    cat("Testing", nrow(sparse_params), "sparse parameter combinations...\n")
    
    results <- vector("list", nrow(sparse_params))
    
    for(i in 1:nrow(sparse_params)) {
        params <- sparse_params[i, ]
        
        # Perform cross-validation for this sparse parameter combination
        cv_result <- .sparse_o2cv(X, Y, nc, nx, ny, params$keepX, params$keepY, 
                                  group, nr_folds)
        
        results[[i]] <- list(
            keepX = params$keepX,
            keepY = params$keepY,
            RMSE = cv_result$RMSE,
            Px = cv_result$Px,
            Py = cv_result$Py,
            Rx = cv_result$Rx,
            Ry = cv_result$Ry
        )
    }
    
    # Find best parameters
    rmse_values <- sapply(results, function(x) x$RMSE)
    best_idx <- which.min(rmse_values)
    best_result <- results[[best_idx]]
    
    cat("Best sparse parameters: keepX =", best_result$keepX, 
        ", keepY =", best_result$keepY, "\n")
    cat("Best RMSE:", best_result$RMSE, "\n")
    
    return(list(
        best_params = best_result,
        all_results = results
    ))
}

#' Internal sparse cross-validation function
#' @keywords internal
.sparse_o2cv <- function(X, Y, nc, nx, ny, keepX, keepY, group, nr_folds) {
    results <- list()
    cls.grp <- getMCCV_cpp(group, n = nr_folds)
    
    fold_results <- vector("list", max(cls.grp, na.rm = TRUE))
    
    for(k in 1:max(cls.grp, na.rm = TRUE)) {
        # Select training group
        Xtr <- X[cls.grp != k & is.na(cls.grp) == FALSE, ]
        Ytr <- Y[cls.grp != k & is.na(cls.grp) == FALSE, ]
        
        # Select testing group
        Xev <- X[cls.grp == k & is.na(cls.grp) == FALSE, ]
        Yev <- Y[cls.grp == k & is.na(cls.grp) == FALSE, ]
        
        tryCatch({
            # Fit sparse O2PLS on training data
            sparse_o2 <- sparse_o2pls(Xtr, Ytr, nc, nx, ny,
                                      keepX = rep(keepX, nc + nx),
                                      keepY = rep(keepY, nc + ny))
            
            # Predict on test data
            Y_pred <- predict(sparse_o2, Xev)
            
            # Calculate errors
            fold_results[[k]] <- list(
                k = k,
                Px = sum((Xev - matrix(0, nrow(Xev), ncol(Xev)))^2),  # Simplified
                Py = sum((Yev - Y_pred)^2),
                Rx = sqrt(mean((Xev - matrix(0, nrow(Xev), ncol(Xev)))^2)),  # Simplified
                Ry = sqrt(mean((Yev - Y_pred)^2))
            )
            
        }, error = function(e) {
            # Handle errors gracefully
            fold_results[[k]] <- list(
                k = k,
                Px = Inf, Py = Inf,
                Rx = Inf, Ry = Inf
            )
        })
    }
    
    # Aggregate results
    valid_results <- fold_results[!sapply(fold_results, is.null)]
    
    if(length(valid_results) == 0) {
        return(list(Px = Inf, Py = Inf, Rx = Inf, Ry = Inf, RMSE = Inf))
    }
    
    Px_mean <- mean(sapply(valid_results, function(x) x$Px), na.rm = TRUE)
    Py_mean <- mean(sapply(valid_results, function(x) x$Py), na.rm = TRUE)
    Rx_mean <- mean(sapply(valid_results, function(x) x$Rx), na.rm = TRUE)
    Ry_mean <- mean(sapply(valid_results, function(x) x$Ry), na.rm = TRUE)
    
    return(list(
        Px = Px_mean,
        Py = Py_mean, 
        Rx = Rx_mean,
        Ry = Ry_mean,
        RMSE = Rx_mean + Ry_mean
    ))
}

#' Stability Selection for Sparse Methods
#' 
#' Performs bootstrap stability selection to identify robust features across
#' multiple bootstrap samples and sparsity levels.
#' 
#' @param X Numeric matrix (samples x variables). The predictor data matrix.
#' @param Y Factor, character vector, or numeric matrix. Response data.
#' @param nc Integer. Number of components to extract.
#' @param keepX_range Vector of integers. Range of keepX values to test (default: seq(10, 100, 10)).
#' @param n_bootstrap Integer. Number of bootstrap samples (default: 100).
#' @param threshold Numeric. Stability threshold between 0 and 1 (default: 0.8).
#'   Variables selected in at least this proportion of bootstrap samples are considered stable.
#' @param method Character. Method to use: "sparse_plsda" or "sparse_o2pls".
#' @param parallel Logical. Use parallel processing (default: TRUE).
#' @param ncores Integer. Number of cores for parallel processing (default: NULL, auto-detect).
#' 
#' @return stability_selection object containing selection probabilities and stable variables.
#' 
#' @details
#' Stability selection is a method for controlling the expected number of false
#' positive selections in high-dimensional variable selection. It works by
#' applying the sparse method to many bootstrap samples and identifying variables
#' that are consistently selected across different samples and sparsity levels.
#' 
#' The threshold parameter controls the stringency of selection. Higher values
#' (closer to 1) result in more conservative selection with fewer false positives
#' but potentially more false negatives.
#' 
#' @examples
#' # Example 1: Stability selection for sparse PLS-DA
#' set.seed(456)
#' n <- 100
#' p <- 150
#' 
#' # Generate data with stable and unstable signals
#' X <- matrix(rnorm(n * p), n, p)
#' 
#' # Create stable signal in first 15 variables
#' stable_signal <- matrix(rnorm(n * 2), n, 2)
#' X[, 1:15] <- stable_signal %*% matrix(rnorm(2 * 15), 2, 15) + 
#'              matrix(rnorm(n * 15, sd = 0.3), n, 15)
#' 
#' # Create unstable/noise signal in variables 16-30
#' X[, 16:30] <- X[, 16:30] + matrix(rnorm(n * 15, sd = 0.8), n, 15)
#' 
#' # Create classes based on stable signal
#' class_probs <- plogis(rowSums(X[, 1:5]) - mean(rowSums(X[, 1:5])))
#' Y <- factor(ifelse(runif(n) < class_probs, "Class1", "Class2"))
#' 
#' colnames(X) <- paste0("Var_", 1:p)
#' 
#' # Perform stability selection (reduced parameters for quick example)
#' stability_result <- stability_selection(
#'     X = X,
#'     Y = Y,
#'     nc = 2,
#'     keepX_range = seq(10, 40, 10),   # Test fewer sparsity levels
#'     n_bootstrap = 50,                # Fewer bootstrap samples for speed
#'     threshold = 0.7,                 # 70% stability threshold
#'     method = "sparse_plsda",
#'     parallel = FALSE,                # Disable parallel for example
#'     ncores = 1
#' )
#' 
#' # View results
#' print(stability_result)
#' 
#' # Extract stable variables for each component
#' stable_comp1 <- stability_result$stable_variables$Component_1
#' stable_comp2 <- stability_result$stable_variables$Component_2
#' 
#' if(!is.null(stable_comp1)) {
#'     cat("Stable variables in Component 1:\n")
#'     print(stable_comp1)
#' }
#' 
#' if(!is.null(stable_comp2)) {
#'     cat("Stable variables in Component 2:\n")
#'     print(stable_comp2)
#' }
#' 
#' # Plot stability results
#' plot(stability_result, type = "heatmap", component = 1)
#' plot(stability_result, type = "barplot", component = 1, top_n = 30)
#' 
#' # Example 2: Stability selection for sparse O2PLS
#' # Generate correlated X and Y matrices
#' Y_matrix <- matrix(rnorm(n * 50), n, 50)
#' # Add correlation with X
#' Y_matrix[, 1:10] <- X[, 1:10] + matrix(rnorm(n * 10, sd = 0.5), n, 10)
#' 
#' colnames(Y_matrix) <- paste0("Y_", 1:50)
#' 
#' # Stability selection for O2PLS
#' \dontrun{
#' stability_o2pls <- stability_selection(
#'     X = X,
#'     Y = Y_matrix,
#'     nc = 2,
#'     keepX_range = seq(15, 45, 15),
#'     n_bootstrap = 100,
#'     threshold = 0.8,
#'     method = "sparse_o2pls",
#'     parallel = TRUE,
#'     ncores = 2
#' )
#' 
#' # Extract stable variables for both X and Y
#' stable_X <- stability_o2pls$stable_variables$Component_1$Variable_Name
#' stable_Y <- stability_o2pls$stable_variables$Component_1$Variable_Name
#' }
#' 
#' # Example 3: Different stability thresholds
#' # Compare conservative vs liberal thresholds
#' \dontrun{
#' # Conservative selection (80% threshold)
#' stability_conservative <- stability_selection(
#'     X, Y, nc = 2, threshold = 0.8, n_bootstrap = 100
#' )
#' 
#' # Liberal selection (60% threshold)  
#' stability_liberal <- stability_selection(
#'     X, Y, nc = 2, threshold = 0.6, n_bootstrap = 100
#' )
#' 
#' # Compare number of stable variables
#' n_stable_cons <- length(stability_conservative$stable_variables$Component_1$Variable)
#' n_stable_lib <- length(stability_liberal$stable_variables$Component_1$Variable)
#' 
#' cat("Conservative (80%):", n_stable_cons, "stable variables\n")
#' cat("Liberal (60%):", n_stable_lib, "stable variables\n")
#' }
#' 
#' @seealso \code{\link{sparse_plsda}}, \code{\link{sparse_o2pls}}, \code{\link{plot.stability_selection}}
#' @author Kai Guo
#' @export
stability_selection <- function(X, Y, nc, keepX_range = seq(10, 100, 10),
                                n_bootstrap = 100, threshold = 0.8,
                                method = "sparse_plsda", parallel = TRUE, 
                                ncores = NULL) {
    
    if(is.null(ncores)) {
        ncores <- min(4, parallel::detectCores() - 1)
    }
    
    n_vars <- ncol(X)
    n_keepX <- length(keepX_range)
    
    # Initialize selection frequency matrix
    selection_freq <- array(0, dim = c(n_vars, n_keepX, nc))
    dimnames(selection_freq) <- list(
        Variables = colnames(X) %||% paste0("Var", 1:n_vars),
        KeepX = paste0("keepX_", keepX_range),
        Components = paste0("Comp", 1:nc)
    )
    
    # Setup parallel processing
    if(parallel && n_bootstrap > 1) {
        cl <- parallel::makeCluster(ncores)
        parallel::clusterEvalQ(cl, library(o2plsda))
        parallel::clusterExport(cl, varlist = c("X", "Y", "nc", "keepX_range", 
                                                "method", "n_vars", "n_keepX"),
                                envir = environment())
    }
    
    # Bootstrap loop
    cat("Performing stability selection with", n_bootstrap, "bootstrap samples...\n")
    pb <- txtProgressBar(min = 0, max = n_bootstrap, style = 3)
    
    bootstrap_results <- if(parallel && n_bootstrap > 1) {
        parallel::parLapply(cl, 1:n_bootstrap, function(b) {
            perform_single_bootstrap(X, Y, nc, keepX_range, method, n_vars, n_keepX)
        })
    } else {
        lapply(1:n_bootstrap, function(b) {
            setTxtProgressBar(pb, b)
            perform_single_bootstrap(X, Y, nc, keepX_range, method, n_vars, n_keepX)
        })
    }
    
    close(pb)
    if(parallel && n_bootstrap > 1) {
        parallel::stopCluster(cl)
    }
    
    # Aggregate bootstrap results
    for(result in bootstrap_results) {
        if(!is.null(result)) {
            selection_freq <- selection_freq + result
        }
    }
    
    # Convert to selection probabilities
    selection_prob <- selection_freq / n_bootstrap
    
    # Identify stable variables
    stable_vars <- identify_stable_variables(selection_prob, threshold, keepX_range)
    
    # Create comprehensive results
    stability_result <- list(
        selection_probabilities = selection_prob,
        stable_variables = stable_vars,
        threshold = threshold,
        n_bootstrap = n_bootstrap,
        keepX_range = keepX_range,
        method = method,
        parameters = list(nc = nc, threshold = threshold)
    )
    
    class(stability_result) <- "stability_selection"
    return(stability_result)
}

# Single bootstrap iteration
perform_single_bootstrap <- function(X, Y, nc, keepX_range, method, n_vars, n_keepX) {
    tryCatch({
        # Bootstrap sample
        n <- nrow(X)
        boot_idx <- sample(1:n, n, replace = TRUE)
        X_boot <- X[boot_idx, , drop = FALSE]
        Y_boot <- Y[boot_idx]
        
        # Initialize selection matrix for this bootstrap
        selection_matrix <- array(0, dim = c(n_vars, n_keepX, nc))
        
        # Test each keepX value
        for(k_idx in 1:n_keepX) {
            keepX_val <- keepX_range[k_idx]
            keepX_vec <- rep(keepX_val, nc)
            
            if(method == "sparse_plsda") {
                fit <- sparse_plsda(X_boot, Y_boot, nc, keepX = keepX_vec, 
                                    tune.keepX = FALSE)
                selected_vars <- fit$selected_vars
            } else if(method == "sparse_o2pls") {
                # Assume Y is matrix for O2PLS
                if(is.factor(Y_boot)) {
                    Y_boot_matrix <- model.matrix(~ . - 1, data = data.frame(class = Y_boot))
                } else {
                    Y_boot_matrix <- as.matrix(Y_boot)
                }
                fit <- sparse_o2pls(X_boot, Y_boot_matrix, nc, nx = 0, ny = 0,
                                    keepX = keepX_vec)
                selected_vars <- fit@sparsity$selected_X
            }
            
            # Record selections for each component
            if(length(selected_vars) > 0) {
                for(comp in 1:nc) {
                    # Determine which variables were selected for this component
                    comp_loadings <- fit$loadings[, comp] # or appropriate loading matrix
                    comp_selected <- which(abs(comp_loadings) > 1e-10)
                    selection_matrix[comp_selected, k_idx, comp] <- 1
                }
            }
        }
        
        return(selection_matrix)
        
    }, error = function(e) {
        return(NULL)  # Skip failed bootstrap samples
    })
}

# Identify stable variables across conditions
identify_stable_variables <- function(selection_prob, threshold, keepX_range) {
    n_vars <- dim(selection_prob)[1]
    n_keepX <- dim(selection_prob)[2] 
    nc <- dim(selection_prob)[3]
    
    stable_results <- list()
    
    for(comp in 1:nc) {
        comp_probs <- selection_prob[, , comp]
        
        # Find variables that are stable across multiple keepX values
        max_prob_per_var <- apply(comp_probs, 1, max)
        stable_vars_comp <- which(max_prob_per_var >= threshold)
        
        if(length(stable_vars_comp) > 0) {
            stable_details <- data.frame(
                Variable = stable_vars_comp,
                Variable_Name = dimnames(selection_prob)[[1]][stable_vars_comp],
                Max_Probability = max_prob_per_var[stable_vars_comp],
                Optimal_KeepX = keepX_range[apply(comp_probs[stable_vars_comp, , drop = FALSE], 
                                                  1, which.max)]
            )
            stable_results[[paste0("Component_", comp)]] <- stable_details
        }
    }
    
    return(stable_results)
}
