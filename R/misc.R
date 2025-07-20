#' @title Summary of an O2PLS object
#'
#' @param object a O2pls object
#' @param ... For consistency
#' @return Detail of O2PLS results
#' @examples 
#' X <- matrix(rnorm(50),10,5)
#' Y <- matrix(rnorm(50),10,5)
#' object <- o2pls(X,Y,1,1,1)
#' summary(object)
#' @author Kai Guo
#' @export
summary.O2pls <- function(object, ...){
    cat("\n######### Summary of the O2PLS results #########\n")
    cat("### Call o2pls(X, Y, nc=",object@params$nc,", nx=",object@params$nx,", ny=",object@params$ny,") ###\n")
    sx <- s2(object@X)
    sy <- s2(object@Y)
    res <- object@results
    d <-data.frame(X=c(res$R2Xcorr,res$R2Xo,1-res$R2X),Y=c(res$R2Ycorr,res$R2Yo,1-res$R2Y))
    varj <- as.data.frame(rbind(res$varXj,res$varYj))
    nc <- object@params$nc
    nx <- object@params$nx
    ny <- object@params$ny
    colnames(varj) <- paste0("LV",1:object@params$nc)
    rownames(varj) <- c("X","Y")
    varx <- rbind(res$varXorth)
    ###if nx =0 or ny =0
    if(nx==0){
        nnx=nc
    }else{
        nnx=nx
    }
    if(ny==0){
        nny=nc
    }else{
        nny=ny
    }
    colnames(varx) <- paste0("LV",1:nnx)
    rownames(varx) <- "X"
    vary <- rbind(res$varYorth)
    colnames(vary) <- paste0("LV",1:nny)
    rownames(vary) <- "Y"
    rownames(d) <- c("Joint","Orthogonal","Noise")
    cat("### Total variation \n")
    cat("### X:",sx, "; Y:",sy," ###\n")
    cat("### Total modeled variation ")
    cat("### X:",round(res$R2X,3), "; Y:", round(res$R2Y,3)," ###\n")
    cat("### Joint, Orthogonal, Noise (proportions) ###\n")
    print(round(d,3))
    cat("### Variation in X joint part predicted by Y Joint part:", round(res$R2Xp,3),"\n")
    cat("### Variation in Y joint part predicted by X Joint part:", round(res$R2Yp,3),"\n")
    cat("### Variation in each Latent Variable (LV) in Joint part: \n")
    print(round(varj,3))
    cat("### Variation in each Latent Variable (LV) in X Orthogonal part: \n")
    print(round(varx,3))
    cat("### Variation in each Latent Variable (LV) in Y Orthogonal part: \n")
    print(round(vary,3))
    cat("\n############################################\n")
}


#' @title Print the summary of O2PLS results.
#' @param x An O2pls object 
#' @param ... For consistency
#' @return NULL
#' @examples 
#' X <- matrix(rnorm(50),10,5)
#' Y <- matrix(rnorm(50),10,5)
#' object <- o2pls(X,Y,1,1,1)
#' print(object)
#' @author Kai Guo
#' @export
print.O2pls <- function (x, ...) {
    summary(x)
}

#' @title Summary of an plsda object
#'
#' @param object a plsda object
#' @param ... For consistency
#' @return Detail of plsda results
#' @examples 
#' X <- matrix(rnorm(500),10,50)
#' Y <- rep(c("a","b"),each=5)
#' fit <- plsda(X,Y,2)
#' summary(fit)
#' @author Kai Guo
#' @export
summary.plsda <- function(object, ...){
    cat("\n######### Summary of the PLS-DA results #########\n")
    cat("### Call plsda(X, Y, nc=",object$nc,") ###\n")
    R2X <- object$xvar
    R2Y <- object$R2Y
    R2Xcum <- cumsum(R2X)
    R2Ycum <- cumsum(R2Y)
    res <- rbind(R2X,R2Xcum,R2Y,R2Ycum)
    print(round(res,3))
    cat("\n############################################\n")
}

#' @title Print the summary of plsda results.
#' @param x An plsda object 
#' @param ... For consistency
#' @return NULL
#' @examples 
#' X <- matrix(rnorm(500),10,50)
#' Y <- rep(c("a","b"),each=5)
#' fit <- plsda(X,Y,2)
#' print(fit)
#' @author Kai Guo
#' @export
print.plsda <- function (x, ...) {
    summary(x)
}

#' Extract the VIP values from the O2PLS-DA object
#' @param x the o2plsda object or plsda object
#' @return a data frame 
#' @export
vip <- function(x){
    res <- x$vip
    return(res)
    
}


#' @title two matrix mutiplication
#' @keywords internal 
#' @return A matrix

eigenmult<-function(A, B){
  return(A %*% B)
}

#' @title three matrix mutiplication
#' @keywords internal 
#' @return A matrix
eigenthree <- function(A, B, C){
    return(A %*% B %*%C)
}


#' @title trans matrix * matrix
#' @keywords internal 
#' @return A matrix
AtA <-function(A) {
    return(crossprod(A))
}

validate_sparse_inputs <- function(X, Y, nc, nx, ny, keepX, keepY) {
    if(!is.matrix(X)) X <- as.matrix(X)
    if(!is.matrix(Y) && !is.factor(Y)) Y <- as.matrix(Y)
    
    if(is.matrix(Y) && nrow(X) != nrow(Y)) {
        stop("X and Y must have the same number of rows")
    }
    if(is.factor(Y) && length(Y) != nrow(X)) {
        stop("Y must have the same length as the number of rows in X")
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
    
    max_components <- min(nrow(X) - 1, ncol(X) - 1)
    if(nc + max(nx, ny) > max_components) {
        stop("Total number of components exceeds data dimensions")
    }
    
    if(!is.null(keepX)) {
        if(!is.numeric(keepX) || any(keepX < 1) || any(keepX > ncol(X))) {
            stop("keepX must contain positive integers not exceeding ncol(X)")
        }
    }
    
    if(!is.null(keepY) && is.matrix(Y)) {
        if(!is.numeric(keepY) || any(keepY < 1) || any(keepY > ncol(Y))) {
            stop("keepY must contain positive integers not exceeding ncol(Y)")
        }
    }
    
    return(TRUE)
}


# Matrix preprocessing
preprocess_matrix <- function(X, scale, center) {
    if(center || scale) {
        return(scale(X, center = center, scale = scale))
    } else {
        return(X)
    }
}



# Calculate various distance metrics
calculate_max_distances <- function(scores, centroids) {
    if(is.null(centroids) || nrow(centroids) == 0) {
        stop("Invalid centroids provided")
    }
    
    distances <- array(NA, dim = c(nrow(scores), nrow(centroids)))
    
    for(i in 1:nrow(centroids)) {
        diff_matrix <- sweep(scores, 2, centroids[i, ], "-")
        distances[, i] <- apply(abs(diff_matrix), 1, max)
    }
    
    return(distances)
}

calculate_centroid_distances <- function(scores, centroids) {
    if(is.null(centroids) || nrow(centroids) == 0) {
        stop("Invalid centroids provided")
    }
    
    distances <- array(NA, dim = c(nrow(scores), nrow(centroids)))
    
    for(i in 1:nrow(centroids)) {
        diff_matrix <- sweep(scores, 2, centroids[i, ], "-")
        distances[, i] <- sqrt(rowSums(diff_matrix^2))
    }
    
    return(distances)
}

calculate_mahalanobis_distances <- function(scores, classification_info) {
    centroids <- classification_info$centroids
    
    if(is.null(centroids) || nrow(centroids) == 0) {
        stop("Invalid centroids provided")
    }
    
    distances <- array(NA, dim = c(nrow(scores), nrow(centroids)))
    
    # Use pooled covariance matrix for Mahalanobis distance
    pooled_cov <- cov(scores)
    
    # Add regularization for numerical stability
    pooled_cov <- pooled_cov + diag(1e-6, ncol(pooled_cov))
    
    tryCatch({
        cov_inv <- solve(pooled_cov)
        
        for(i in 1:nrow(centroids)) {
            diff_matrix <- sweep(scores, 2, centroids[i, ], "-")
            distances[, i] <- sqrt(diag(diff_matrix %*% cov_inv %*% t(diff_matrix)))
        }
    }, error = function(e) {
        # Fall back to Euclidean distance
        warning("Mahalanobis distance calculation failed, using Euclidean distance")
        distances <- calculate_centroid_distances(scores, centroids)
    })
    
    return(distances)
}

# Calculate prediction probabilities from distances
calculate_prediction_probabilities <- function(distances) {
    # Convert distances to probabilities using softmax
    # Smaller distances should give higher probabilities
    neg_distances <- -distances
    
    # Numerical stability: subtract max from each row
    max_vals <- apply(neg_distances, 1, max, na.rm = TRUE)
    neg_distances_stable <- neg_distances - max_vals
    
    exp_distances <- exp(neg_distances_stable)
    probabilities <- exp_distances / rowSums(exp_distances, na.rm = TRUE)
    
    # Handle any NaN values
    probabilities[is.nan(probabilities)] <- 1 / ncol(probabilities)
    
    return(probabilities)
}

calculate_classification_error <- function(true_labels, predicted_labels, measure) {
    if(length(true_labels) != length(predicted_labels)) {
        return(1.0)
    }
    
    if(measure == "BER") {
        confusion <- table(true_labels, predicted_labels)
        if(nrow(confusion) == 0 || ncol(confusion) == 0) return(1.0)
        
        class_errors <- numeric(nrow(confusion))
        for(i in 1:nrow(confusion)) {
            if(sum(confusion[i, ]) > 0) {
                class_errors[i] <- 1 - confusion[i, i] / sum(confusion[i, ])
            } else {
                class_errors[i] <- 1.0
            }
        }
        return(mean(class_errors, na.rm = TRUE))
    } else if(measure == "accuracy") {
        return(mean(true_labels == predicted_labels, na.rm = TRUE))
    } else {
        return(mean(true_labels != predicted_labels, na.rm = TRUE))
    }
}

calculate_regression_error <- function(true_values, predicted_values, measure) {
    if(length(true_values) != length(predicted_values)) {
        return(Inf)
    }
    
    if(is.matrix(true_values)) true_values <- as.vector(true_values)
    if(is.matrix(predicted_values)) predicted_values <- as.vector(predicted_values)
    
    if(measure == "RMSE") {
        return(sqrt(mean((true_values - predicted_values)^2, na.rm = TRUE)))
    } else if(measure == "R2") {
        ss_res <- sum((true_values - predicted_values)^2, na.rm = TRUE)
        ss_tot <- sum((true_values - mean(true_values, na.rm = TRUE))^2, na.rm = TRUE)
        if(ss_tot == 0) return(0)
        return(max(0, 1 - ss_res/ss_tot))
    } else if(measure == "MAE") {
        return(mean(abs(true_values - predicted_values), na.rm = TRUE))
    } else {
        return(sqrt(mean((true_values - predicted_values)^2, na.rm = TRUE)))
    }
}

#' Summary method for sparse O2PLS results
#' @export
summary.SparseO2pls <- function(object, ...) {
    # Handle both S4 and list objects
    if(inherits(object, "SparseO2pls")) {
        call_obj <- object@call
        params <- object@params
        results <- object@results
        sparsity <- object@sparsity
        X_data <- object@X
        Y_data <- object@Y
    } else {
        call_obj <- object$call
        params <- object$params
        results <- object$results
        sparsity <- object$sparsity
        X_data <- object$X
        Y_data <- object$Y
    }
    
    cat("\n")
    cat("######### Summary of Sparse O2PLS Results #########\n")
    cat("### Call: ")
    if(!is.null(call_obj)) {
        cat(deparse(call_obj), "\n")
    } else {
        cat("sparse_o2pls(...)\n")
    }
    
    # Parameters
    cat("### Parameters ###\n")
    cat("   Components: nc =", params$nc)
    if(!is.null(params$nx) && params$nx > 0) cat(", nx =", params$nx)
    if(!is.null(params$ny) && params$ny > 0) cat(", ny =", params$ny)
    cat("\n")
    
    keepX_str <- if(!is.null(params$keepX)) paste(params$keepX, collapse = ", ") else "auto"
    keepY_str <- if(!is.null(params$keepY)) paste(params$keepY, collapse = ", ") else "auto"
    
    cat("   Sparsity: keepX = [", keepX_str, "]")
    cat(", keepY = [", keepY_str, "]\n")
    
    if(!is.null(params$penalty)) {
        penalty_str <- params$penalty
        lambda_x_str <- if(!is.null(params$lambda_x)) params$lambda_x else "default"
        lambda_y_str <- if(!is.null(params$lambda_y)) params$lambda_y else "default"
        cat("   Penalty:", penalty_str, "(λ_X =", lambda_x_str, 
            ", λ_Y =", lambda_y_str, ")\n")
    }
    
    # Data dimensions
    cat("### Data Dimensions ###\n")
    cat("   X (samples × variables):", nrow(X_data), "×", ncol(X_data), "\n")
    cat("   Y (samples × variables):", nrow(Y_data), "×", ncol(Y_data), "\n")
    
    # Variable selection results
    cat("### Variable Selection ###\n")
    n_selected_X <- length(sparsity$selected_X)
    n_selected_Y <- length(sparsity$selected_Y)
    total_X <- ncol(X_data)
    total_Y <- ncol(Y_data)
    
    cat("   Selected X variables:", n_selected_X, "/", total_X, 
        " (", round(100 * n_selected_X / total_X, 1), "% kept)\n")
    cat("   Selected Y variables:", n_selected_Y, "/", total_Y, 
        " (", round(100 * n_selected_Y / total_Y, 1), "% kept)\n")
    
    sparsity_X <- if(!is.null(sparsity$sparsity_X)) sparsity$sparsity_X else (1 - n_selected_X/total_X)
    sparsity_Y <- if(!is.null(sparsity$sparsity_Y)) sparsity$sparsity_Y else (1 - n_selected_Y/total_Y)
    
    cat("   Sparsity achieved: X =", round(sparsity_X * 100, 1), "%",
        ", Y =", round(sparsity_Y * 100, 1), "%\n")
    
    # Variance explained
    cat("### Variance Explained (Sparse Model) ###\n")
    
    R2X <- if(!is.null(results$R2X)) results$R2X else NA
    R2Y <- if(!is.null(results$R2Y)) results$R2Y else NA
    R2X_joint <- if(!is.null(results$R2X_joint)) results$R2X_joint else NA
    R2Y_joint <- if(!is.null(results$R2Y_joint)) results$R2Y_joint else NA
    
    if(!is.na(R2X) && !is.na(R2Y)) {
        cat("   Total variance: X =", round(R2X, 3), ", Y =", round(R2Y, 3), "\n")
    }
    
    if(!is.na(R2X_joint) && !is.na(R2Y_joint)) {
        cat("   Joint variance: X =", round(R2X_joint, 3), ", Y =", round(R2Y_joint, 3), "\n")
    }
    
    # Orthogonal variance (if components exist)
    if((!is.null(params$nx) && params$nx > 0) || (!is.null(params$ny) && params$ny > 0)) {
        R2X_orth <- if(!is.null(results$R2X_orth)) results$R2X_orth else 0
        R2Y_orth <- if(!is.null(results$R2Y_orth)) results$R2Y_orth else 0
        cat("   Orthogonal variance: X =", round(R2X_orth, 3),
            ", Y =", round(R2Y_orth, 3), "\n")
    }
    
    # Predictability
    if(params$nc > 0) {
        R2X_pred <- if(!is.null(results$R2X_pred)) results$R2X_pred else NA
        R2Y_pred <- if(!is.null(results$R2Y_pred)) results$R2Y_pred else NA
        
        if(!is.na(R2X_pred) && !is.na(R2Y_pred)) {
            cat("   Predictability: X =", round(R2X_pred, 3),
                ", Y =", round(R2Y_pred, 3), "\n")
        }
    }
    
    cat("###################################################\n")
    cat("\n")
    
    invisible(list(
        parameters = params,
        dimensions = list(X = dim(X_data), Y = dim(Y_data)),
        selection = list(
            selected_X = n_selected_X, total_X = total_X,
            selected_Y = n_selected_Y, total_Y = total_Y,
            sparsity_X = sparsity_X,
            sparsity_Y = sparsity_Y
        ),
        variance = list(
            R2X = R2X, R2Y = R2Y,
            R2X_joint = R2X_joint, R2Y_joint = R2Y_joint,
            R2X_orth = if(exists("R2X_orth")) R2X_orth else NULL,
            R2Y_orth = if(exists("R2Y_orth")) R2Y_orth else NULL
        )
    ))
}

#' Print method for sparse O2PLS results
#' @export
print.SparseO2pls <- function(x, ...) {
    summary(x, ...)
    invisible(x)
}


#' @title Summary of sparse PLS-DA results  
#' @export
summary.sparse_plsda <- function(object, ...) {
    cat("\n######### Summary of Sparse PLS-DA Results #########\n")
    cat("### Components:", object@ncomp, "###\n")
    cat("### KeepX:", paste(object@keepX, collapse = ", "), "###\n")
    cat("### Selected variables:", length(object@selected_vars), 
        "out of", nrow(object@loadings), "###\n")
    cat("### Classes:", paste(object@Y_levels, collapse = ", "), "###\n")
    
    if(!is.null(object@classification$overall_error)) {
        cat("### Classification error:", 
            round(object@classification$overall_error * 100, 2), "% ###\n")
    }
    
    cat("\n####################################################\n")
}

#' @export
print.sparse_plsda <- function(x, ...) {
    summary.sparse_plsda(x)
    invisible(x)
}
#' @export
print.sparse_o2pls <- function(x, ...) {
    if(inherits(x, "SparseO2pls")) {
        summary.SparseO2pls(x, ...)
    } else {
        summary.sparse_o2pls(x, ...)
    }
    invisible(x)
}

#' @export
summary.sparse_o2pls <- function(object, ...) {
    # Handle both S4 and list objects
    if(inherits(object, "SparseO2pls")) {
        call_obj <- object@call
        params <- object@params
        sparsity <- object@sparsity
        X_data <- object@X
        Y_data <- object@Y
    } else {
        call_obj <- object$call
        params <- object$params
        sparsity <- object$sparsity
        X_data <- object$X
        Y_data <- object$Y
    }
    
    cat("\n######### Summary of Sparse O2PLS Results #########\n")
    cat("### Components: nc =", params$nc)
    if(!is.null(params$nx) && params$nx > 0) cat(", nx =", params$nx)
    if(!is.null(params$ny) && params$ny > 0) cat(", ny =", params$ny)
    cat("\n")
    
    n_selected_X <- length(sparsity$selected_X)
    n_selected_Y <- length(sparsity$selected_Y)
    total_X <- ncol(X_data)
    total_Y <- ncol(Y_data)
    
    cat("### Selected X variables:", n_selected_X, "/", total_X, "\n")
    cat("### Selected Y variables:", n_selected_Y, "/", total_Y, "\n")
    cat("###################################################\n")
}
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
create_stratified_folds <- function(Y, n_folds) {
    if(is.factor(Y)) {
        # Enhanced stratified sampling for classification
        folds <- integer(length(Y))
        for(level in levels(Y)) {
            idx <- which(Y == level)
            n_level <- length(idx)
            if(n_level >= n_folds) {
                folds[idx] <- sample(rep(1:n_folds, length.out = n_level))
            } else {
                # For small classes, assign to different folds
                folds[idx] <- sample(1:n_folds, n_level, replace = FALSE)
            }
        }
        return(folds)
    } else {
        # Simple random sampling for regression
        return(sample(rep(1:n_folds, length.out = length(Y))))
    }
}
find_optimal_parameters <- function(results, measure) {
    # Extract scores and filter out failed results
    scores <- sapply(results, function(x) {
        if(is.null(x$mean_score) || !is.finite(x$mean_score)) {
            return(if(measure %in% c("BER", "RMSE", "MAE")) Inf else -Inf)
        }
        return(x$mean_score)
    })
    
    # Find optimal parameters
    if(measure %in% c("BER", "RMSE", "MAE")) {
        # Lower is better
        optimal_idx <- which.min(scores)
    } else {
        # Higher is better (accuracy, R2)
        optimal_idx <- which.max(scores)
    }
    
    if(length(optimal_idx) == 0) {
        stop("No valid results found during tuning")
    }
    
    return(results[[optimal_idx]])
}


`%||%` <- function(x, y) if(is.null(x)) y else x
