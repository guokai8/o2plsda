#' @title Partial least squares discriminant analysis
#' @description Perform a PLS discriminant analysis
#' @importFrom stats model.matrix cor 
#' @param X a matrix of predictor variables.
#' @param Y a single vector indicate the group
#' @param nc the number of pls components
#' @param scale logical indicating whether X must be scaled (suggest TRUE).
#' @param center logical indicating whether X must be centered (suggest TRUE).
#' @param cv logical indicating whether cross-validation will be performed or not (suggest TRUE).
#' @param nr_folds Integer to indicate the folds for cross validation.
#' @param tol Convergence tolerance (default: 1e-6)
#' @param max_iter Maximum number of iterations (default: 100)
#' @param q2_threshold Q2 threshold for component selection (default: 0.05)
#' @return a list containing PLS-DA results
#' @export

plsda <- function(X, Y, nc, scale = TRUE, center = TRUE, cv = TRUE, 
                  nr_folds = 5, tol = 1e-6, max_iter = 100, q2_threshold = 0.05) {
    # Input validation
    X <- as.matrix(X)
    if (!is.numeric(nc) || nc < 1 || nc > min(nrow(X), ncol(X))) {
        stop("nc must be a positive integer less than min(nrow(X), ncol(X))")
    }
    if (!is.numeric(nr_folds) || nr_folds < 2) {
        stop("nr_folds must be >= 2")
    }
    
    # Convert Y to model matrix (dummy coding)
    Y <- model.matrix(~ -1 + Y)
    
    if (nrow(X) != nrow(Y)) {
        stop("Y should have same length as number of X rows")
    }
    
    # Scale X (Y scaling is optional for discriminant analysis)
    if (isTRUE(scale)) {
        X <- scale(X, center, scale = TRUE)
    } else if (isTRUE(center)) {
        X <- scale(X, center, scale = FALSE)
    }
    
    # Fit PLS model
    fit <- .pls(X, Y, nc, cv = cv, nr_folds = nr_folds, 
                tol = tol, max_iter = max_iter, q2_threshold = q2_threshold)
    
    # Process results
    nc <- fit$nc
    score <- fit$Tt
    colnames(score) <- paste0("LV", 1:nc)
    rownames(score) <- rownames(X)
    
    Xloading <- fit$W
    colnames(Xloading) <- paste0("LV", 1:nc)
    rownames(Xloading) <- colnames(X)
    
    # Calculate VIP scores more efficiently
    W <- Xloading
    p <- ncol(X)
    Th <- score
    VIP <- matrix(0, nrow = p, ncol = nc)
    cor2 <- cor(Y, Th, use = "pairwise.complete.obs")^2
    
    # Vectorized VIP calculation
    for (h in 1:nc) {
        R <- colSums(cor2[, 1:h, drop = FALSE])
        VIP[, h] <- sqrt(p * (R %*% t(W[, 1:h, drop = FALSE]^2))/sum(R))
    }
    
    rownames(VIP) <- rownames(W)
    colnames(VIP) <- paste0("comp", 1:nc)
    
    # Prepare return values
    xvar <- fit$x_var
    names(xvar) <- paste0('LV', 1:nc)
    R2Y <- fit$R2Y
    names(R2Y) <- paste0('LV', 1:nc)
    
    res <- list(
        nc = nc,
        scores = score,
        Xloading = Xloading,
        vip = VIP,
        xvar = xvar,
        R2Y = R2Y,
        PRESS = fit$PRESS,
        Q2 = fit$Q2,
        Q2G = fit$Q2G
    )
    
    class(res) <- "plsda"
    return(res)
}

#' @title Internal PLS function
#' @description Internal function to perform PLS calculations
#' @keywords internal
.pls <- function(X, Y, nc, cv = TRUE, nr_folds = 5, 
                 tol = 1e-6, max_iter = 100, q2_threshold = 0.05) {
    n <- nrow(X)
    m <- ncol(X)
    q <- ncol(Y)
    px <- ncol(X)
    Xt <- X
    Yt <- Y
    xvar <- sum(Xt^2)
    x_var <- NULL
    
    # Determine optimal number of components if cv=TRUE
    if (isTRUE(cv)) {
        Xs <- svd(Xt, nu = 0, nv = 0)
        rank_X <- sum(Xs$d > tol)
        nc <- min(nc, n, rank_X)
    }
    
    # Initialize matrices
    W <- matrix(0, m, nc)
    U <- matrix(0, n, nc)
    Tt <- matrix(0, n, nc)
    Ch <- matrix(0, q, nc)
    Ph <- matrix(0, px, nc)
    RSS <- rbind(rep(n-1, q), matrix(NA, nc, q))
    PRESS <- matrix(NA, nc, q)
    Q2 <- matrix(NA, nc, q)
    
    # Create cross-validation folds
    fold <- split(sample(1:n), rep(1:nr_folds, length = n))
    
    # PLS2 algorithm with improved convergence checks
    for (i in 1:nc) {
        u <- Yt[, 1]
        wo <- rep(1, m)
        
        # Inner loop for convergence
        for (iter in 1:max_iter) {
            w <- t(Xt) %*% u / sum(u^2)
            w <- w / sqrt(sum(w^2))
            tp <- Xt %*% w
            ct <- t(Yt) %*% tp / sum(tp^2)
            u <- Yt %*% ct / sum(ct^2)
            
            if (sum((w - wo)^2) < tol) break
            wo <- w
        }
        
        p <- t(Xt) %*% tp / sum(tp^2)
        
        # Cross-validation
        RSS[i+1, ] <- colSums((Yt - tp %*% t(ct))^2)
        
        if (isTRUE(cv)) {
            press <- matrix(0, n, q)
            
            for (k in 1:nr_folds) {
                omit <- fold[[k]]
                uh.cv <- Yt[-omit, 1]
                wh.cvt <- rep(1, px)
                
                # Cross-validation inner loop
                for (itcv in 1:max_iter) {
                    wh.cv <- t(Xt[-omit, ]) %*% uh.cv / sum(uh.cv^2)
                    wh.cv <- wh.cv / sqrt(sum(wh.cv^2))
                    th.cv <- Xt[-omit, ] %*% wh.cv
                    ch.cv <- t(Yt[-omit, ]) %*% th.cv / sum(th.cv^2)
                    uh.cv <- Yt[-omit, ] %*% ch.cv / sum(ch.cv^2)
                    
                    if (sum((wh.cv - wh.cvt)^2) < tol) break
                    wh.cvt <- wh.cv
                }
                
                Yhat.cv <- (Xt[omit, ] %*% wh.cv) %*% t(ch.cv)
                press[omit, ] <- (Yt[omit, ] - Yhat.cv)^2
            }
            
            PRESS[i, ] <- colSums(press)
            Q2[i, ] <- 1 - PRESS[i, ]/RSS[i, ]
        }
        
        # Deflation
        Xt <- Xt - (tp %*% t(p))
        Yt <- Yt - (tp %*% t(ct))
        
        # Store results
        W[, i] <- w
        U[, i] <- u
        Tt[, i] <- tp
        Ch[, i] <- ct
        Ph[, i] <- p
    }
    
    # Calculate Q2 global and determine optimal number of components
    Q2G <- if (cv) 1 - rowSums(PRESS)/rowSums(RSS[-nc, ]) else NULL
    
    if (isTRUE(cv)) {
        ns <- max(which(Q2G >= q2_threshold))
        if (ns < nc) {
            nc <- max(2, ns)  # Ensure at least 2 components
            # Truncate matrices to optimal number of components
            W <- W[, 1:nc, drop = FALSE]
            U <- U[, 1:nc, drop = FALSE]
            Ph <- Ph[, 1:nc, drop = FALSE]
            Tt <- Tt[, 1:nc, drop = FALSE]
            Ch <- Ch[, 1:nc, drop = FALSE]
        }
    }
    
    # Calculate explained variance
    x_var <- sapply(1:nc, function(j) {
        v <- t(Tt[, j, drop = FALSE]) %*% X
        as.numeric(v %*% t(v) / crossprod(Tt[, j], Tt[, j]) / xvar)
    })
    
    # Calculate R2Y
    R2Y <- colMeans(cor(Y, Tt)^2)
    
    return(list(
        nc = nc,
        W = W,
        U = U,
        Tt = Tt,
        Ch = Ch,
        Ph = Ph,
        x_var = x_var,
        R2Y = R2Y,
        PRESS = PRESS,
        Q2 = Q2,
        Q2G = Q2G
    ))
}
