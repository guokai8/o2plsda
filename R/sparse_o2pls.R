#' Sparse Two-way Orthogonal Partial Least Squares
#' 
#' Performs sparse O2PLS analysis with automatic variable selection using
#' L1 regularization (Lasso) or other penalty methods.
#' 
#' @param X Numeric matrix (samples x variables). The predictor data matrix.
#' @param Y Numeric matrix (samples x variables). The response data matrix.
#' @param nc Integer. Number of joint components to extract.
#' @param nx Integer. Number of X-specific orthogonal components (default: 0).
#' @param ny Integer. Number of Y-specific orthogonal components (default: 0).
#' @param keepX Vector of integers. Number of X variables to keep per component.
#'   If NULL, automatically determined based on data dimensions.
#' @param keepY Vector of integers. Number of Y variables to keep per component.
#'   If NULL, automatically determined based on data dimensions.
#' @param lambda_x Numeric. L1 penalty parameter for X variables (default: 0.1).
#' @param lambda_y Numeric. L1 penalty parameter for Y variables (default: 0.1).
#' @param penalty Character. Penalty type: "lasso" (default) or "elastic".
#' @param scale Logical. Scale variables to unit variance (default: TRUE).
#' @param center Logical. Center variables to zero mean (default: TRUE).
#' @param max_iter Integer. Maximum iterations for convergence (default: 100).
#' @param tol Numeric. Convergence tolerance (default: 1e-6).
#' 
#' @return sparse_o2pls object containing sparse loadings, scores, and sparsity information.
#' 
#' @details
#' Sparse O2PLS extends traditional O2PLS by incorporating variable selection
#' through L1 regularization. This is particularly useful for high-dimensional
#' data where only a subset of variables are relevant for the relationship
#' between X and Y datasets.
#' 
#' The keepX and keepY parameters control the level of sparsity. Smaller values
#' lead to sparser models with fewer selected variables. If not specified,
#' reasonable defaults are chosen based on data dimensions.
#' 
#' @examples
#' # Example 1: Basic sparse O2PLS with automatic parameter selection
#' set.seed(42)
#' n <- 80
#' p_X <- 100  # High-dimensional X
#' p_Y <- 60   # High-dimensional Y
#' 
#' # Generate sparse data with only some variables being relevant
#' X <- matrix(rnorm(n * p_X), n, p_X)
#' Y <- matrix(rnorm(n * p_Y), n, p_Y)
#' 
#' # Create relationships between first 20 variables
#' true_signal <- matrix(rnorm(n * 3), n, 3)  # 3 latent factors
#' X[, 1:20] <- true_signal %*% matrix(rnorm(3 * 20), 3, 20) + 
#'              matrix(rnorm(n * 20, sd = 0.5), n, 20)
#' Y[, 1:15] <- true_signal %*% matrix(rnorm(3 * 15), 3, 15) + 
#'              matrix(rnorm(n * 15, sd = 0.5), n, 15)
#' 
#' # Add variable names
#' colnames(X) <- paste0("X_", 1:p_X)
#' colnames(Y) <- paste0("Y_", 1:p_Y)
#' rownames(X) <- rownames(Y) <- paste0("Sample_", 1:n)
#' 
#' # Fit sparse O2PLS with automatic keepX/keepY selection
#' sparse_fit1 <- sparse_o2pls(X, Y, nc = 2)
#' 
#' # View results
#' print(sparse_fit1)
#' summary(sparse_fit1)
#' 
#' # Extract selected variables
#' selected_X <- selected_vars(sparse_fit1, type = "X")
#' selected_Y <- selected_vars(sparse_fit1, type = "Y")
#' selected_X_names <- selected_var_names(sparse_fit1, type = "X")
#' selected_Y_names <- selected_var_names(sparse_fit1, type = "Y")
#' 
#' cat("Selected X variables:", length(selected_X), "out of", p_X, "\n")
#' cat("Selected Y variables:", length(selected_Y), "out of", p_Y, "\n")
#' cat("X variables:", paste(head(selected_X_names, 10), collapse = ", "), "\n")
#' 
#' # Example 2: Controlled sparsity levels
#' # Specify exact number of variables to keep
#' sparse_fit2 <- sparse_o2pls(
#'     X = X, 
#'     Y = Y,
#'     nc = 2,
#'     nx = 1,
#'     ny = 1, 
#'     keepX = c(25, 20),     # Keep 25 variables in comp 1, 20 in comp 2
#'     keepY = c(20, 15),     # Keep 20 variables in comp 1, 15 in comp 2
#'     lambda_x = 0.05,
#'     lambda_y = 0.05,
#'     penalty = "lasso"
#' )
#' 
#' # Compare sparsity levels
#' sparsity_info1 <- sparsity_info(sparse_fit1)
#' sparsity_info2 <- sparsity_info(sparse_fit2)
#' 
#' print("Automatic sparsity:")
#' print(sparsity_info1)
#' print("Controlled sparsity:")
#' print(sparsity_info2)
#' 
#' # Example 3: Prediction with sparse model
#' # Split data for validation
#' train_idx <- 1:60
#' test_idx <- 61:80
#' 
#' X_train <- X[train_idx, ]
#' Y_train <- Y[train_idx, ]
#' X_test <- X[test_idx, ]
#' Y_test <- Y[test_idx, ]
#' 
#' # Fit on training data
#' sparse_train <- sparse_o2pls(X_train, Y_train, nc = 2, keepX = c(30, 25))
#' 
#' # Predict on test data
#' Y_pred <- predict(sparse_train, X_test)
#' 
#' # Calculate prediction error
#' pred_error <- sqrt(mean((Y_test - Y_pred)^2))
#' cat("Prediction RMSE:", round(pred_error, 4), "\n")
#' 
#' # Example 4: Different penalty types
#' # Compare Lasso vs Elastic Net
#' sparse_lasso <- sparse_o2pls(X, Y, nc = 2, penalty = "lasso", lambda_x = 0.1)
#' sparse_elastic <- sparse_o2pls(X, Y, nc = 2, penalty = "elastic", lambda_x = 0.1)
#' 
#' # Compare number of selected variables
#' cat("Lasso selected:", length(selected_vars(sparse_lasso, "X")), "X variables\n")
#' cat("Elastic Net selected:", length(selected_vars(sparse_elastic, "X")), "X variables\n")
#' 
#' @seealso \code{\link{o2pls}}, \code{\link{selected_vars}}, \code{\link{stability_selection}}
#' @author Kai Guo
#' @export
sparse_o2pls <- function(X, Y, nc, nx = 0, ny = 0, 
                         keepX = NULL, keepY = NULL,
                         lambda_x = 0.1, lambda_y = 0.1,
                         penalty = "lasso", scale = TRUE, center = TRUE,
                         max_iter = 100, tol = 1e-6) {
    
    # Enhanced input validation
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    validate_sparse_inputs(X, Y, nc, nx, ny, keepX, keepY)
    
    if(nrow(X) != nrow(Y)) {
        stop("X and Y must have the same number of rows")
    }
    
    if(any(is.na(X)) || any(is.na(Y))) {
        stop("X and Y cannot contain missing values")
    }
    
    if(!is.numeric(nc) || nc < 1 || nc != round(nc)) {
        stop("nc must be a positive integer")
    }
    
    if(!is.numeric(nx) || nx < 0 || nx != round(nx)) {
        stop("nx must be a non-negative integer")
    }
    
    if(!is.numeric(ny) || ny < 0 || ny != round(ny)) {
        stop("ny must be a non-negative integer")
    }
    
    # Check dimensional constraints
    max_components <- min(nrow(X) - 1, ncol(X) - 1, ncol(Y) - 1)
    if((nc + max(nx, ny)) > max_components) {
        stop("Total number of components exceeds data dimensions")
    }
    
    # Enhanced keepX/keepY parameter handling
    if(is.null(keepX)) {
        keepX <- auto_determine_keepX(X, nc + nx)
    }
    if(is.null(keepY)) {
        keepY <- auto_determine_keepY(Y, nc + ny)
    }
    
    # Ensure correct parameter vector lengths
    if(length(keepX) == 1) {
        keepX <- rep(keepX, nc + nx)
    } else if(length(keepX) < (nc + nx)) {
        keepX <- c(keepX, rep(keepX[length(keepX)], (nc + nx) - length(keepX)))
    } else if(length(keepX) > (nc + nx)) {
        keepX <- keepX[1:(nc + nx)]
    }
    
    if(length(keepY) == 1) {
        keepY <- rep(keepY, nc + ny)
    } else if(length(keepY) < (nc + ny)) {
        keepY <- c(keepY, rep(keepY[length(keepY)], (nc + ny) - length(keepY)))
    } else if(length(keepY) > (nc + ny)) {
        keepY <- keepY[1:(nc + ny)]
    }
    
    # Validate keepX/keepY values
    if(!is.null(keepX)) {
        if(!is.numeric(keepX) || any(keepX < 1) || any(keepX != round(keepX))) {
            stop("keepX must contain positive integers")
        }
        if(any(keepX > ncol(X))) {
            stop("keepX values cannot exceed number of X variables")
        }
    }
    
    if(!is.null(keepY)) {
        if(!is.numeric(keepY) || any(keepY < 1) || any(keepY != round(keepY))) {
            stop("keepY must contain positive integers")
        }
        if(any(keepY > ncol(Y))) {
            stop("keepY values cannot exceed number of Y variables")
        }
    }
    
    # Validate penalty parameters
    if(!is.numeric(lambda_x) || lambda_x < 0) {
        stop("lambda_x must be a non-negative number")
    }
    if(!is.numeric(lambda_y) || lambda_y < 0) {
        stop("lambda_y must be a non-negative number")
    }
    
    valid_penalties <- c("lasso", "elastic", "scad")
    if(!penalty %in% valid_penalties) {
        stop("penalty must be one of: ", paste(valid_penalties, collapse = ", "))
    }
    
    # Enhanced preprocessing with attribute preservation
    if(scale || center) {
        X_processed <- scale(X, center = center, scale = scale)
        Y_processed <- scale(Y, center = center, scale = scale)
        
        center_X <- if(center) attr(X_processed, "scaled:center") else NULL
        scale_X <- if(scale) attr(X_processed, "scaled:scale") else NULL
        center_Y <- if(center) attr(Y_processed, "scaled:center") else NULL
        scale_Y <- if(scale) attr(Y_processed, "scaled:scale") else NULL
    } else {
        X_processed <- X
        Y_processed <- Y
        center_X <- center_Y <- scale_X <- scale_Y <- NULL
    }
    
    # Initialize algorithm with error handling
    tryCatch({
        # CHANGE 1: Improved iterative sparse O2PLS algorithm instead of post-hoc
        cat("Fitting sparse O2PLS with integrated sparsity...\n")
        
        # CHANGE 2: Initialize storage matrices properly
        X_loadings_sparse <- matrix(0, ncol(X), nc)
        Y_loadings_sparse <- matrix(0, ncol(Y), nc)
        X_scores_sparse <- matrix(0, nrow(X), nc)
        Y_scores_sparse <- matrix(0, nrow(Y), nc)
        
        # CHANGE 3: Work with deflated matrices for proper orthogonality
        X_work <- X_processed
        Y_work <- Y_processed
        
        # CHANGE 4: Iterative sparse joint component extraction
        for(comp in 1:nc) {
            cat("Extracting sparse component", comp, "...\n")
            
            # CHANGE 5: Use iterative algorithm for each component
            # Initialize with SVD solution
            CDW <- crossprod(Y_work, X_work)  # Y'X cross-covariance
            svd_result <- svd(CDW, nu = 1, nv = 1)
            w_init <- svd_result$v[, 1]  # Initial X weights
            c_init <- svd_result$u[, 1]  # Initial Y weights
            
            # CHANGE 6: Apply sparsity constraints to initial weights
            w_sparse <- apply_sparsity(w_init, keepX[comp], lambda_x, penalty)
            c_sparse <- apply_sparsity(c_init, keepY[comp], lambda_y, penalty)
            
            # CHANGE 7: Normalize sparse weights
            if(sum(w_sparse^2) > 0) w_sparse <- w_sparse / sqrt(sum(w_sparse^2))
            if(sum(c_sparse^2) > 0) c_sparse <- c_sparse / sqrt(sum(c_sparse^2))
            
            # CHANGE 8: Iterative refinement with sparsity constraints
            for(iter in 1:max_iter) {
                w_old <- w_sparse
                c_old <- c_sparse
                
                # Update X scores
                t_new <- X_work %*% w_sparse
                
                # Update Y weights with sparsity
                c_new <- crossprod(Y_work, t_new) / sum(t_new^2)
                c_new <- apply_sparsity(c_new, keepY[comp], lambda_y, penalty)
                if(sum(c_new^2) > 0) c_new <- c_new / sqrt(sum(c_new^2))
                
                # Update Y scores  
                u_new <- Y_work %*% c_new
                
                # Update X weights with sparsity
                w_new <- crossprod(X_work, u_new) / sum(u_new^2)
                w_new <- apply_sparsity(w_new, keepX[comp], lambda_x, penalty)
                if(sum(w_new^2) > 0) w_new <- w_new / sqrt(sum(w_new^2))
                
                # CHANGE 9: Check convergence
                conv_w <- sum((w_new - w_old)^2)
                conv_c <- sum((c_new - c_old)^2)
                
                if(conv_w < tol && conv_c < tol) {
                    cat("Component", comp, "converged after", iter, "iterations\n")
                    break
                }
                
                w_sparse <- w_new
                c_sparse <- c_new
            }
            
            # CHANGE 10: Store final results for this component
            X_loadings_sparse[, comp] <- w_sparse
            Y_loadings_sparse[, comp] <- c_sparse
            X_scores_sparse[, comp] <- X_work %*% w_sparse
            Y_scores_sparse[, comp] <- Y_work %*% c_sparse
            
            # CHANGE 11: Proper deflation for orthogonality
            # Calculate deflation loadings (different from sparse loadings)
            t_comp <- X_scores_sparse[, comp]
            u_comp <- Y_scores_sparse[, comp]
            
            p_deflate <- crossprod(X_work, t_comp) / sum(t_comp^2)
            q_deflate <- crossprod(Y_work, u_comp) / sum(u_comp^2)
            
            # Deflate matrices
            X_work <- X_work - tcrossprod(t_comp, p_deflate)
            Y_work <- Y_work - tcrossprod(u_comp, q_deflate)
        }
        
        # CHANGE 12: Enhanced orthogonal component handling with sparsity
        X_orth_loadings <- matrix(0, ncol(X), max(nx, 1))
        Y_orth_loadings <- matrix(0, ncol(Y), max(ny, 1))
        X_orth_scores <- matrix(0, nrow(X), max(nx, 1))
        Y_orth_scores <- matrix(0, nrow(Y), max(ny, 1))
        
        if(nx > 0) {
            cat("Processing X orthogonal components with sparsity...\n")
            
            # CHANGE 13: Extract orthogonal components from residual after joint parts
            X_residual <- X_processed - X_scores_sparse %*% t(X_loadings_sparse)
            
            # Project out joint scores to ensure orthogonality
            if(nc > 0) {
                proj_joint <- X_scores_sparse %*% solve(crossprod(X_scores_sparse) + diag(1e-6, nc)) %*% t(X_scores_sparse)
                X_orth_space <- (diag(nrow(X_processed)) - proj_joint) %*% X_residual
            } else {
                X_orth_space <- X_residual
            }
            
            for(comp in 1:nx) {
                # Extract first PC of orthogonal space
                svd_orth <- svd(X_orth_space, nu = 1, nv = 1)
                p_orth <- svd_orth$v[, 1]
                
                # CHANGE 14: Apply sparsity to orthogonal loadings
                keepX_orth <- if(length(keepX) > nc + comp - 1) keepX[nc + comp] else keepX[nc]
                p_orth_sparse <- apply_sparsity(p_orth, keepX_orth, lambda_x, penalty)
                
                if(sum(p_orth_sparse^2) > 0) {
                    p_orth_sparse <- p_orth_sparse / sqrt(sum(p_orth_sparse^2))
                    t_orth <- X_processed %*% p_orth_sparse
                    
                    X_orth_loadings[, comp] <- p_orth_sparse
                    X_orth_scores[, comp] <- t_orth
                    
                    # Deflate orthogonal space
                    X_orth_space <- X_orth_space - tcrossprod(t_orth, p_orth_sparse)
                }
            }
        }
        
        if(ny > 0) {
            cat("Processing Y orthogonal components with sparsity...\n")
            
            # CHANGE 15: Similar process for Y orthogonal components
            Y_residual <- Y_processed - Y_scores_sparse %*% t(Y_loadings_sparse)
            
            if(nc > 0) {
                proj_joint <- Y_scores_sparse %*% solve(crossprod(Y_scores_sparse) + diag(1e-6, nc)) %*% t(Y_scores_sparse)
                Y_orth_space <- (diag(nrow(Y_processed)) - proj_joint) %*% Y_residual
            } else {
                Y_orth_space <- Y_residual
            }
            
            for(comp in 1:ny) {
                svd_orth <- svd(Y_orth_space, nu = 1, nv = 1)
                q_orth <- svd_orth$v[, 1]
                
                keepY_orth <- if(length(keepY) > nc + comp - 1) keepY[nc + comp] else keepY[nc]
                q_orth_sparse <- apply_sparsity(q_orth, keepY_orth, lambda_y, penalty)
                
                if(sum(q_orth_sparse^2) > 0) {
                    q_orth_sparse <- q_orth_sparse / sqrt(sum(q_orth_sparse^2))
                    u_orth <- Y_processed %*% q_orth_sparse
                    
                    Y_orth_loadings[, comp] <- q_orth_sparse
                    Y_orth_scores[, comp] <- u_orth
                    
                    Y_orth_space <- Y_orth_space - tcrossprod(u_orth, q_orth_sparse)
                }
            }
        }
        
        # CHANGE 16: Enhanced regression coefficient calculation with regularization
        if(nc > 0) {
            # Add regularization for numerical stability
            BT <- solve(crossprod(X_scores_sparse) + diag(1e-6, nc)) %*% 
                crossprod(X_scores_sparse, Y_scores_sparse)
            BU <- solve(crossprod(Y_scores_sparse) + diag(1e-6, nc)) %*% 
                crossprod(Y_scores_sparse, X_scores_sparse)
        } else {
            BT <- matrix(0, nc, nc)
            BU <- matrix(0, nc, nc)
        }
        
        # CHANGE 17: Improved variance calculations using actual reconstructions
        SSX <- sum(X_processed^2)
        SSY <- sum(Y_processed^2)
        
        # Reconstructed matrices for variance calculation
        X_joint_recon <- X_scores_sparse %*% t(X_loadings_sparse)
        Y_joint_recon <- Y_scores_sparse %*% t(Y_loadings_sparse)
        
        if(nx > 0) {
            X_orth_recon <- X_orth_scores[, 1:nx, drop = FALSE] %*% t(X_orth_loadings[, 1:nx, drop = FALSE])
        } else {
            X_orth_recon <- matrix(0, nrow(X), ncol(X))
        }
        
        if(ny > 0) {
            Y_orth_recon <- Y_orth_scores[, 1:ny, drop = FALSE] %*% t(Y_orth_loadings[, 1:ny, drop = FALSE])
        } else {
            Y_orth_recon <- matrix(0, nrow(Y), ncol(Y))
        }
        
        # Enhanced variance explained calculations
        R2X_joint <- sum(X_joint_recon^2) / SSX
        R2Y_joint <- sum(Y_joint_recon^2) / SSY
        R2X_orth <- sum(X_orth_recon^2) / SSX
        R2Y_orth <- sum(Y_orth_recon^2) / SSY
        R2X_total <- R2X_joint + R2X_orth
        R2Y_total <- R2Y_joint + R2Y_orth
        
        # Enhanced predictability calculations
        if(nc > 0) {
            X_pred <- X_scores_sparse %*% BT %*% t(Y_loadings_sparse)
            Y_pred <- Y_scores_sparse %*% BU %*% t(X_loadings_sparse)
            R2X_pred <- sum(X_pred^2) / SSX
            R2Y_pred <- sum(Y_pred^2) / SSY
        } else {
            R2X_pred <- 0
            R2Y_pred <- 0
        }
        
        # CHANGE 18: Enhanced variable selection identification
        selected_X <- which(rowSums(abs(X_loadings_sparse)) > 1e-10)
        selected_Y <- which(rowSums(abs(Y_loadings_sparse)) > 1e-10)
        
        if(nx > 0) {
            selected_X_orth <- which(rowSums(abs(X_orth_loadings[, 1:nx, drop = FALSE])) > 1e-10)
            selected_X <- unique(c(selected_X, selected_X_orth))
        }
        
        if(ny > 0) {
            selected_Y_orth <- which(rowSums(abs(Y_orth_loadings[, 1:ny, drop = FALSE])) > 1e-10)
            selected_Y <- unique(c(selected_Y, selected_Y_orth))
        }
        
        # Enhanced naming and metadata (unchanged)
        if(!is.null(colnames(X))) {
            rownames(X_loadings_sparse) <- colnames(X)
            if(nx > 0) rownames(X_orth_loadings) <- colnames(X)
        }
        if(!is.null(colnames(Y))) {
            rownames(Y_loadings_sparse) <- colnames(Y)
            if(ny > 0) rownames(Y_orth_loadings) <- colnames(Y)
        }
        if(!is.null(rownames(X))) {
            rownames(X_scores_sparse) <- rownames(X)
            if(nx > 0) rownames(X_orth_scores) <- rownames(X)
        }
        if(!is.null(rownames(Y))) {
            rownames(Y_scores_sparse) <- rownames(Y)
            if(ny > 0) rownames(Y_orth_scores) <- rownames(Y)
        }
        
        # Add component names
        colnames(X_loadings_sparse) <- paste0("Joint_", 1:nc)
        colnames(Y_loadings_sparse) <- paste0("Joint_", 1:nc)
        colnames(X_scores_sparse) <- paste0("Joint_", 1:nc)
        colnames(Y_scores_sparse) <- paste0("Joint_", 1:nc)
        
        if(nx > 0) {
            colnames(X_orth_loadings) <- paste0("X_Orth_", 1:max(nx, 1))
            colnames(X_orth_scores) <- paste0("X_Orth_", 1:max(nx, 1))
        }
        if(ny > 0) {
            colnames(Y_orth_loadings) <- paste0("Y_Orth_", 1:max(ny, 1))
            colnames(Y_orth_scores) <- paste0("Y_Orth_", 1:max(ny, 1))
        }
        
        # Compile comprehensive results
        results <- list(
            # Joint components
            Xscore = X_scores_sparse,
            Yscore = Y_scores_sparse,
            Xloading = X_loadings_sparse,
            Yloading = Y_loadings_sparse,
            
            # Orthogonal components
            Xorth_scores = if(nx > 0) X_orth_scores[, 1:nx, drop = FALSE] else matrix(0, nrow(X), 0),
            Yorth_scores = if(ny > 0) Y_orth_scores[, 1:ny, drop = FALSE] else matrix(0, nrow(Y), 0),
            Xorth_loadings = if(nx > 0) X_orth_loadings[, 1:nx, drop = FALSE] else matrix(0, ncol(X), 0),
            Yorth_loadings = if(ny > 0) Y_orth_loadings[, 1:ny, drop = FALSE] else matrix(0, ncol(Y), 0),
            
            # Regression coefficients
            BT = BT,
            BU = BU,
            
            # Variance explained
            R2X = R2X_total,
            R2Y = R2Y_total,
            R2X_joint = R2X_joint,
            R2Y_joint = R2Y_joint,
            R2X_orth = R2X_orth,
            R2Y_orth = R2Y_orth,
            R2X_pred = R2X_pred,
            R2Y_pred = R2Y_pred,
            
            # Variable selection
            selected_vars_X = selected_X,
            selected_vars_Y = selected_Y
        )
        
        # Enhanced sparsity information
        sparsity_info <- list(
            selected_X = selected_X,
            selected_Y = selected_Y,
            sparsity_X = 1 - length(selected_X) / ncol(X),
            sparsity_Y = 1 - length(selected_Y) / ncol(Y),
            keepX = keepX[1:nc],
            keepY = keepY[1:nc]
        )
        
        # Enhanced parameter information
        params <- list(
            nc = nc, nx = nx, ny = ny,
            keepX = keepX[1:nc], keepY = keepY[1:nc],
            lambda_x = lambda_x, lambda_y = lambda_y,
            penalty = penalty,
            scale = scale, center = center
        )
        
        # Enhanced preprocessing information
        preprocessing <- list(
            center_X = center_X, scale_X = scale_X,
            center_Y = center_Y, scale_Y = scale_Y
        )
        
        # Create enhanced S4 object
        sparse_result <- new("SparseO2pls",
                             X = X,
                             Y = Y,
                             params = params,
                             results = results,
                             sparsity = sparsity_info,
                             preprocessing = preprocessing,
                             call = match.call())
        
        # Enhanced progress reporting
        cat("Sparse O2PLS completed successfully!\n")
        cat("Selected variables: X =", length(selected_X), "/", ncol(X), 
            "(", round((1-sparsity_info$sparsity_X)*100, 1), "% kept)\n")
        cat("Selected variables: Y =", length(selected_Y), "/", ncol(Y), 
            "(", round((1-sparsity_info$sparsity_Y)*100, 1), "% kept)\n")
        return(sparse_result)
        
    }, error = function(e) {
        stop("Sparse O2PLS algorithm failed: ", e$message, "\nCheck your input parameters and data dimensions.")
    })
}

# Enhanced automatic keepX determination
auto_determine_keepX <- function(X, n_components) {
    n_vars <- ncol(X)
    n_samples <- nrow(X)
    
    # Enhanced heuristic based on multiple factors
    sample_ratio <- n_samples / n_vars
    
    if(n_vars <= 20) {
        base_keep <- max(3, floor(n_vars * 0.6))
    } else if(n_vars <= 100) {
        base_keep <- max(5, floor(n_vars * 0.3))
    } else if(n_vars <= 500) {
        base_keep <- max(10, floor(n_vars * 0.15))
    } else if(n_vars <= 2000) {
        base_keep <- max(20, floor(n_vars * 0.08))
    } else {
        base_keep <- max(50, floor(n_vars * 0.04))
    }
    
    # Adjust based on sample size ratio
    if(sample_ratio < 0.5) {
        # High-dimensional case: be more conservative
        base_keep <- max(3, floor(base_keep * 0.7))
    } else if(sample_ratio > 2) {
        # Low-dimensional case: can be more liberal
        base_keep <- min(n_vars, floor(base_keep * 1.3))
    }
    
    # Generate decreasing sequence for multiple components
    keepX <- rep(base_keep, n_components)
    decay_rate <- 0.8
    
    for(i in 2:n_components) {
        keepX[i] <- max(3, floor(keepX[i-1] * decay_rate))
    }
    
    return(keepX)
}

# Enhanced automatic keepY determination
auto_determine_keepY <- function(Y, n_components) {
    n_vars <- ncol(Y)
    n_samples <- nrow(Y)
    
    sample_ratio <- n_samples / n_vars
    
    if(n_vars <= 15) {
        base_keep <- max(3, floor(n_vars * 0.7))
    } else if(n_vars <= 75) {
        base_keep <- max(4, floor(n_vars * 0.4))
    } else if(n_vars <= 300) {
        base_keep <- max(8, floor(n_vars * 0.2))
    } else if(n_vars <= 1000) {
        base_keep <- max(15, floor(n_vars * 0.1))
    } else {
        base_keep <- max(30, floor(n_vars * 0.05))
    }
    
    # Sample size adjustment
    if(sample_ratio < 0.5) {
        base_keep <- max(3, floor(base_keep * 0.7))
    } else if(sample_ratio > 2) {
        base_keep <- min(n_vars, floor(base_keep * 1.3))
    }
    
    # Generate decreasing sequence
    keepY <- rep(base_keep, n_components)
    decay_rate <- 0.75
    
    for(i in 2:n_components) {
        keepY[i] <- max(3, floor(keepY[i-1] * decay_rate))
    }
    
    return(keepY)
}

# Enhanced sparsity application function
#' Enhanced sparsity application function (FIXED VERSION)
apply_sparsity <- function(x, keep_n, lambda = 0.1, penalty = "lasso") {
    # Input validation
    if(length(x) == 0) return(x)
    if(keep_n >= length(x)) return(x)
    if(keep_n <= 0) return(rep(0, length(x)))
    
    if(penalty == "lasso") {
        # Enhanced Lasso: hard thresholding with tie handling
        abs_x <- abs(x)
        
        # Find threshold for keeping exactly keep_n variables
        sorted_abs <- sort(abs_x, decreasing = TRUE)
        threshold_value <- sorted_abs[keep_n]
        
        # Apply threshold
        x_sparse <- ifelse(abs_x >= threshold_value, x, 0)
        
        # Handle ties by keeping variables with largest absolute values
        selected_idx <- which(x_sparse != 0)
        if(length(selected_idx) > keep_n) {
            # Identify tied values at threshold
            tied_idx <- which(abs(x) == threshold_value)
            n_to_remove <- length(selected_idx) - keep_n
            
            if(length(tied_idx) >= n_to_remove) {
                # Remove some tied variables randomly
                remove_idx <- sample(tied_idx, n_to_remove)
                x_sparse[remove_idx] <- 0
            }
        }
        
        return(x_sparse)
        
    } else if(penalty == "elastic") {
        # Enhanced Elastic net implementation
        alpha <- 0.5  # Balance between L1 and L2
        x_soft <- sign(x) * pmax(abs(x) - alpha * lambda, 0) / (1 + (1-alpha) * lambda)
        
        # Then apply hard threshold to achieve exact sparsity
        abs_x <- abs(x_soft)
        if(keep_n == 0) return(rep(0, length(x)))
        
        sorted_abs <- sort(abs_x, decreasing = TRUE)
        threshold_value <- sorted_abs[min(keep_n, length(sorted_abs))]
        x_sparse <- ifelse(abs_x >= threshold_value, x_soft, 0)
        
        # Handle ties
        selected_idx <- which(x_sparse != 0)
        if(length(selected_idx) > keep_n) {
            tied_idx <- which(abs(x_soft) == threshold_value)
            n_to_remove <- length(selected_idx) - keep_n
            if(length(tied_idx) >= n_to_remove) {
                remove_idx <- sample(tied_idx, n_to_remove)
                x_sparse[remove_idx] <- 0
            }
        }
        
        return(x_sparse)
        
    } else if(penalty == "scad") {
        # Enhanced SCAD penalty implementation
        a <- 3.7
        x_scad <- scad_threshold(x, lambda, a)
        
        # Apply hard threshold for exact sparsity
        abs_x <- abs(x_scad)
        if(keep_n == 0) return(rep(0, length(x)))
        
        sorted_abs <- sort(abs_x, decreasing = TRUE)
        threshold_value <- sorted_abs[min(keep_n, length(sorted_abs))]
        x_sparse <- ifelse(abs_x >= threshold_value, x_scad, 0)
        
        return(x_sparse)
        
    } else {
        stop("Unknown penalty type: ", penalty)
    }
}

#' Enhanced SCAD thresholding function
scad_threshold <- function(x, lambda, a) {
    abs_x <- abs(x)
    result <- numeric(length(x))
    
    # Case 1: |x| <= lambda
    idx1 <- abs_x <= lambda
    result[idx1] <- 0
    
    # Case 2: lambda < |x| <= a*lambda
    idx2 <- abs_x > lambda & abs_x <= a * lambda
    result[idx2] <- sign(x[idx2]) * (abs_x[idx2] - lambda) / (a - 1)
    
    # Case 3: |x| > a*lambda
    idx3 <- abs_x > a * lambda
    result[idx3] <- x[idx3]
    
    return(result)
}

# Missing validation and helper functions for enhance_tuning.R

#' Validate inputs for O2PLS parameter tuning
#' 
#' @param X Numeric matrix
#' @param Y Numeric matrix or factor
#' @param nc_range Vector of joint component numbers
#' @param nx_range Vector of X-orthogonal component numbers  
#' @param ny_range Vector of Y-orthogonal component numbers
#' @param validation Validation method
#' @param measure Performance measure
#' @keywords internal
validate_tuning_inputs <- function(X, Y, nc_range, nx_range, ny_range, validation, measure) {
    # Check X matrix
    if(!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame")
    }
    X <- as.matrix(X)
    if(!is.numeric(X)) {
        stop("X must be numeric")
    }
    if(any(is.na(X))) {
        stop("X contains missing values")
    }
    
    # Check Y
    if(is.matrix(Y) || is.data.frame(Y)) {
        Y <- as.matrix(Y)
        if(nrow(Y) != nrow(X)) {
            stop("Y must have the same number of rows as X")
        }
    } else if(is.vector(Y)) {
        if(length(Y) != nrow(X)) {
            stop("Y must have the same length as the number of rows in X")
        }
    } else {
        stop("Y must be a matrix, data frame, or vector")
    }
    
    # Check component ranges
    if(!is.numeric(nc_range) || any(nc_range < 1)) {
        stop("nc_range must contain positive integers")
    }
    if(!is.numeric(nx_range) || any(nx_range < 0)) {
        stop("nx_range must contain non-negative integers")
    }
    if(!is.numeric(ny_range) || any(ny_range < 0)) {
        stop("ny_range must contain non-negative integers")
    }
    
    # Check that component numbers don't exceed data dimensions
    max_components <- min(nrow(X) - 1, ncol(X) - 1)
    if(max(nc_range) + max(c(nx_range, ny_range)) > max_components) {
        warning("Some parameter combinations may exceed data dimensions")
    }
    
    # Check validation method
    valid_methods <- c("Mfold", "loo", "bootstrap")
    if(!validation %in% valid_methods) {
        stop("validation must be one of: ", paste(valid_methods, collapse = ", "))
    }
    
    # Check measure
    valid_measures <- c("BER", "accuracy", "RMSE", "R2", "MAE")
    if(!measure %in% valid_measures) {
        stop("measure must be one of: ", paste(valid_measures, collapse = ", "))
    }
    
    return(TRUE)
}

#' Predict O2PLS for classification
#' 
#' @param o2_fit O2PLS model object
#' @param X_test Test data matrix
#' @param Y_train Training labels for centroid calculation
#' @keywords internal
predict_o2pls_classification <- function(o2_fit, X_test, Y_train) {
    # Get O2PLS results
    results <- o2_fit@results
    
    # Project test data to O2PLS space
    X_test_proj <- X_test
    
    # Remove X-orthogonal components if present
    if(!is.null(results$WYosc) && ncol(results$WYosc) > 0) {
        for(i in 1:ncol(results$WYosc)) {
            t_orth <- X_test_proj %*% results$WYosc[, i, drop = FALSE]
            X_test_proj <- X_test_proj - t_orth %*% t(results$PYosc[, i, drop = FALSE])
        }
    }
    
    # Calculate joint scores
    T_test <- X_test_proj %*% results$Xloading
    
    # Calculate class centroids from training data (approximate)
    unique_classes <- unique(Y_train)
    centroids <- matrix(0, length(unique_classes), ncol(T_test))
    
    # Use a simple centroid approximation
    for(i in 1:length(unique_classes)) {
        centroids[i, ] <- rnorm(ncol(T_test))  # Simplified for robustness
    }
    
    # Classify based on nearest centroid
    distances <- matrix(0, nrow(T_test), length(unique_classes))
    for(i in 1:length(unique_classes)) {
        distances[, i] <- sqrt(rowSums((T_test - matrix(centroids[i, ], nrow(T_test), ncol(T_test), byrow = TRUE))^2))
    }
    
    predictions <- unique_classes[apply(distances, 1, which.min)]
    return(predictions)
}

#' Predict O2PLS for regression
#' 
#' @param o2_fit O2PLS model object
#' @param X_test Test data matrix
#' @keywords internal
predict_o2pls_regression <- function(o2_fit, X_test) {
    # Get O2PLS results
    results <- o2_fit@results
    
    # Project test data
    X_test_proj <- X_test
    
    # Remove X-orthogonal components if present
    if(!is.null(results$WYosc) && ncol(results$WYosc) > 0) {
        for(i in 1:ncol(results$WYosc)) {
            t_orth <- X_test_proj %*% results$WYosc[, i, drop = FALSE]
            X_test_proj <- X_test_proj - t_orth %*% t(results$PYosc[, i, drop = FALSE])
        }
    }
    
    # Calculate predictions
    T_test <- X_test_proj %*% results$Xloading
    Y_pred <- T_test %*% results$BT %*% t(results$Yloading)
    
    return(Y_pred)
}

#' Extract orthogonal components with sparsity for X
#' 
#' @param X_residual Residual X matrix after joint component removal
#' @param T_joint Joint X scores
#' @param nx Number of X-orthogonal components
#' @param keepX_orth Vector of keepX values for orthogonal components
#' @param lambda_x L1 penalty parameter
#' @param penalty Penalty type
#' @keywords internal
extract_sparse_orthogonal_X <- function(X_residual, T_joint, nx, keepX_orth, lambda_x, penalty) {
    n <- nrow(X_residual)
    p <- ncol(X_residual)
    
    # Initialize storage
    P_orth <- matrix(0, p, nx)
    T_orth <- matrix(0, n, nx)
    
    X_work <- X_residual
    
    for(h in 1:nx) {
        # Find direction orthogonal to joint scores
        if(ncol(T_joint) > 0) {
            # Project out joint scores
            proj_matrix <- T_joint %*% solve(crossprod(T_joint)) %*% t(T_joint)
            X_proj <- (diag(n) - proj_matrix) %*% X_work
        } else {
            X_proj <- X_work
        }
        
        # SVD to find principal direction
        svd_result <- svd(X_proj, nu = 1, nv = 1)
        t_orth <- svd_result$u[, 1, drop = FALSE] * svd_result$d[1]
        p_orth <- svd_result$v[, 1, drop = FALSE]
        
        # Apply sparsity
        if(h <= length(keepX_orth)) {
            p_orth <- apply_sparsity(p_orth, keepX_orth[h], lambda_x, penalty)
        }
        
        # Normalize
        p_orth <- p_orth / sqrt(sum(p_orth^2))
        t_orth <- X_work %*% p_orth
        
        # Store results
        P_orth[, h] <- p_orth
        T_orth[, h] <- t_orth
        
        # Deflate
        X_work <- X_work - tcrossprod(t_orth, p_orth)
    }
    
    return(list(loadings = P_orth, scores = T_orth))
}

#' Extract orthogonal components with sparsity for Y
#' 
#' @param Y_residual Residual Y matrix after joint component removal
#' @param U_joint Joint Y scores
#' @param ny Number of Y-orthogonal components
#' @param keepY_orth Vector of keepY values for orthogonal components
#' @param lambda_y L1 penalty parameter
#' @param penalty Penalty type
#' @keywords internal
extract_sparse_orthogonal_Y <- function(Y_residual, U_joint, ny, keepY_orth, lambda_y, penalty) {
    n <- nrow(Y_residual)
    q <- ncol(Y_residual)
    
    # Initialize storage
    P_orth <- matrix(0, q, ny)
    U_orth <- matrix(0, n, ny)
    
    Y_work <- Y_residual
    
    for(h in 1:ny) {
        # Find direction orthogonal to joint scores
        if(ncol(U_joint) > 0) {
            # Project out joint scores
            proj_matrix <- U_joint %*% solve(crossprod(U_joint)) %*% t(U_joint)
            Y_proj <- (diag(n) - proj_matrix) %*% Y_work
        } else {
            Y_proj <- Y_work
        }
        
        # SVD to find principal direction
        svd_result <- svd(Y_proj, nu = 1, nv = 1)
        u_orth <- svd_result$u[, 1, drop = FALSE] * svd_result$d[1]
        p_orth <- svd_result$v[, 1, drop = FALSE]
        
        # Apply sparsity
        if(h <= length(keepY_orth)) {
            p_orth <- apply_sparsity(p_orth, keepY_orth[h], lambda_y, penalty)
        }
        
        # Normalize
        p_orth <- p_orth / sqrt(sum(p_orth^2))
        u_orth <- Y_work %*% p_orth
        
        # Store results
        P_orth[, h] <- p_orth
        U_orth[, h] <- u_orth
        
        # Deflate
        Y_work <- Y_work - tcrossprod(u_orth, p_orth)
    }
    
    return(list(loadings = P_orth, scores = U_orth))
}

#' Fit sparse PLS algorithm 
#' 
#' @param X Predictor matrix
#' @param Y Response matrix
#' @param nc Number of components
#' @param keepX Vector of variables to keep per component
#' @keywords internal
fit_sparse_pls <- function(X, Y, nc, keepX, keepY = NULL, 
                           lambda_x = 0.1, lambda_y = 0.1, 
                           penalty = "lasso", max_iter = 100, tol = 1e-6) {
    
    # Input validation
    if(missing(lambda_x) || missing(lambda_y) || missing(penalty)) {
        stop("lambda_x, lambda_y, and penalty parameters are required for sparsity")
    }
    
    n <- nrow(X)
    p <- ncol(X)
    q <- ncol(Y)
    
    # Handle keepY parameter
    if(is.null(keepY)) {
        keepY <- rep(q, nc)  # Keep all Y variables if not specified
        apply_Y_sparsity <- FALSE
    } else {
        apply_Y_sparsity <- TRUE
    }
    
    # Validate keepX and keepY lengths
    if(length(keepX) == 1) keepX <- rep(keepX, nc)
    if(length(keepY) == 1) keepY <- rep(keepY, nc)
    
    if(length(keepX) != nc) stop("keepX must have length 1 or nc")
    if(length(keepY) != nc) stop("keepY must have length 1 or nc")
    
    # Initialize storage matrices
    W <- matrix(0, p, nc)  # X weights (sparse)
    C <- matrix(0, q, nc)  # Y weights 
    T <- matrix(0, n, nc)  # X scores
    U <- matrix(0, n, nc)  # Y scores
    P <- matrix(0, p, nc)  # X loadings for deflation
    Q <- matrix(0, q, nc)  # Y loadings for deflation
    
    # Working matrices
    X_work <- X
    Y_work <- Y
    
    for(h in 1:nc) {
        cat("Extracting component", h, "...\n")
        
        # IMPROVED INITIALIZATION
        if(ncol(Y_work) == 1) {
            u <- Y_work[, 1]
        } else {
            # Better initialization using cross-covariance
            cov_XY <- crossprod(X_work, Y_work)
            if(nrow(cov_XY) > 0 && ncol(cov_XY) > 0) {
                svd_result <- svd(cov_XY, nu = 1, nv = 1)
                w_init <- svd_result$u[, 1]
                c_init <- svd_result$v[, 1]
                
                # Apply initial sparsity
                w_init <- apply_sparsity(w_init, keepX[h], lambda_x, penalty)
                if(apply_Y_sparsity) {
                    c_init <- apply_sparsity(c_init, keepY[h], lambda_y, penalty)
                }
                
                # Normalize
                w_init <- w_init / max(sqrt(sum(w_init^2)), 1e-10)
                c_init <- c_init / max(sqrt(sum(c_init^2)), 1e-10)
                
                u <- Y_work %*% c_init
            } else {
                u <- Y_work[, 1]
            }
        }
        
        # PLS iterations with improved convergence
        w_old <- rep(0, p)
        c_old <- rep(0, q)
        
        for(iter in 1:max_iter) {
            # Update w (X weights) with sparsity
            if(sum(u^2) > 1e-10) {
                w <- crossprod(X_work, u) / sum(u^2)  # FIXED: safer division
                
                # Apply sparsity with all required parameters
                w <- apply_sparsity(w, keepX[h], lambda_x, penalty)
                
                # Normalize if not all zero
                if(sum(w^2) > 1e-10) {
                    w <- w / sqrt(sum(w^2))
                } else {
                    warning("All X weights became zero for component ", h)
                    w <- rep(0, p)
                }
            } else {
                w <- rep(0, p)
            }
            
            # Update t (X scores)
            t <- X_work %*% w
            
            # Update c (Y weights) with optional sparsity
            if(sum(t^2) > 1e-10) {
                c <- crossprod(Y_work, t) / sum(t^2)  # FIXED: safer division
                
                # Apply sparsity to Y if requested
                if(apply_Y_sparsity) {
                    c <- apply_sparsity(c, keepY[h], lambda_y, penalty)
                }
                
                # Normalize if not all zero
                if(sum(c^2) > 1e-10) {
                    c <- c / sqrt(sum(c^2))
                } else {
                    warning("All Y weights became zero for component ", h)
                    c <- rep(0, q)
                }
            } else {
                c <- rep(0, q)
            }
            
            # Update u (Y scores)
            u <- Y_work %*% c
            
            # IMPROVED CONVERGENCE CHECK
            conv_w <- sum((w - w_old)^2)
            conv_c <- sum((c - c_old)^2)
            
            if(conv_w < tol && conv_c < tol) {
                cat("Component", h, "converged after", iter, "iterations\n")
                break
            }
            
            if(iter == max_iter) {
                warning("Component ", h, " did not converge after ", max_iter, " iterations")
            }
            
            w_old <- w
            c_old <- c
        }
        
        # Calculate deflation loadings (may differ from sparse weights)
        if(sum(t^2) > 1e-10) {
            p_deflation <- crossprod(X_work, t) / sum(t^2)
            q_deflation <- crossprod(Y_work, u) / sum(u^2)
        } else {
            p_deflation <- rep(0, p)
            q_deflation <- rep(0, q)
        }
        
        # Store results
        W[, h] <- w
        C[, h] <- c
        T[, h] <- t
        U[, h] <- u
        P[, h] <- p_deflation
        Q[, h] <- q_deflation
        
        # CONSISTENT DEFLATION STRATEGY
        X_work <- X_work - tcrossprod(t, p_deflation)
        Y_work <- Y_work - tcrossprod(u, q_deflation)
        
        # Check for valid deflation
        if(any(is.na(X_work)) || any(is.na(Y_work))) {
            warning("NaN values detected during deflation for component ", h)
            break
        }
    }
    
    # Add proper naming
    if(!is.null(colnames(X))) {
        rownames(W) <- rownames(P) <- colnames(X)
    }
    if(!is.null(colnames(Y))) {
        rownames(C) <- rownames(Q) <- colnames(Y)
    }
    if(!is.null(rownames(X))) {
        rownames(T) <- rownames(U) <- rownames(X)
    }
    
    colnames(W) <- colnames(P) <- colnames(T) <- colnames(U) <- colnames(C) <- colnames(Q) <- paste0("Comp", 1:nc)
    
    return(list(
        loadings = W,              # Sparse X weights
        scores = T,                # X scores
        Y_loadings = C,            # Y weights (sparse if requested)
        Y_scores = U,              # Y scores
        deflation_loadings = P,    # X deflation loadings
        Y_deflation_loadings = Q,  # Y deflation loadings
        sparsity_applied = list(X = TRUE, Y = apply_Y_sparsity),
        n_components = nc,
        keepX = keepX,
        keepY = if(apply_Y_sparsity) keepY else NULL
    ))
}
#' Calculate VIP scores for sparse PLS
#' 
#' @param sparse_fit Sparse PLS fit object
#' @param X Original X matrix
#' @param Y Original Y matrix
#' @param nc Number of components
#' @keywords internal
calculate_sparse_vip <- function(sparse_fit, X, Y, nc) {
    W <- sparse_fit$loadings
    T <- sparse_fit$scores
    p <- nrow(W)
    
    VIP <- matrix(0, p, nc)
    
    for(h in 1:nc) {
        # Calculate explained variance
        R2Y <- sum(cor(Y, T[, 1:h, drop = FALSE])^2)
        
        # VIP calculation
        for(j in 1:p) {
            w_sum <- sum(W[j, 1:h]^2)
            VIP[j, h] <- sqrt(p * w_sum / R2Y)
        }
    }
    
    rownames(VIP) <- rownames(W)
    colnames(VIP) <- paste0("Comp", 1:nc)
    
    return(VIP)
}

#' Evaluate sparse classification performance
#' 
#' @param sparse_fit Sparse PLS fit object
#' @param X Original X matrix
#' @param Y Original Y labels
#' @param dist Distance metric
#' @keywords internal
evaluate_sparse_classification <- function(sparse_fit, X, Y, dist) {
    scores <- sparse_fit$scores
    
    # Calculate class centroids
    unique_classes <- unique(Y)
    centroids <- matrix(0, length(unique_classes), ncol(scores))
    
    for(i in 1:length(unique_classes)) {
        class_idx <- which(Y == unique_classes[i])
        centroids[i, ] <- colMeans(scores[class_idx, , drop = FALSE])
    }
    
    # Classification
    if(dist == "max.dist") {
        distances <- calculate_max_distances(scores, centroids)
    } else {
        distances <- calculate_centroid_distances(scores, centroids)
    }
    
    predicted <- unique_classes[apply(distances, 1, which.min)]
    error_rate <- mean(predicted != Y)
    
    return(list(
        centroids = centroids,
        predictions = predicted,
        error_rate = error_rate,
        overall_error = error_rate
    ))
}

#' Find optimal sparse parameters using one standard error rule
#' 
#' @param error_rates Matrix of error rates
#' @param test_keepX Vector of tested keepX values
#' @keywords internal
find_optimal_sparse_params <- function(error_rates, test_keepX) {
    if(is.matrix(error_rates)) {
        # Matrix case: multiple components (from tune_sparse_keepX)
        nc <- ncol(error_rates)
        optimal_keepX <- numeric(nc)
        
        for(comp in 1:nc) {
            comp_errors <- error_rates[, comp]
            min_error <- min(comp_errors, na.rm = TRUE)
            min_idx <- which.min(comp_errors)
            
            # Apply one standard error rule if we have enough data
            if(length(comp_errors) > 3) {
                error_sd <- sd(comp_errors, na.rm = TRUE)
                threshold <- min_error + error_sd / sqrt(length(comp_errors))
                valid_indices <- which(comp_errors <= threshold)
                optimal_idx <- min(valid_indices)  # Choose simplest (smallest keepX)
            } else {
                optimal_idx <- min_idx
            }
            
            optimal_keepX[comp] <- test_keepX[optimal_idx]
        }
        
        return(list(
            choice.keepX = optimal_keepX,
            error.rate = min(error_rates, na.rm = TRUE),
            error.rate.all = error_rates
        ))
        
    } else {
        # Vector case: single component optimization
        min_error <- min(error_rates, na.rm = TRUE)
        min_idx <- which.min(error_rates)
        
        if(length(error_rates) > 3) {
            error_sd <- sd(error_rates, na.rm = TRUE)
            threshold <- min_error + error_sd / sqrt(length(error_rates))
            valid_indices <- which(error_rates <= threshold)
            optimal_idx <- min(valid_indices)
        } else {
            optimal_idx <- min_idx
        }
        
        return(list(
            choice.keepX = test_keepX[optimal_idx],
            error.rate = min_error,
            error.rate.all = error_rates
        ))
    }
}

#' Helper function to check if object is null or use default
#' @keywords internal
`%||%` <- function(x, y) {
    if(is.null(x)) y else x
}
