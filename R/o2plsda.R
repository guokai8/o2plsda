#' @title Orthogonal partial least squares discriminant analysis
#' @description Computes orthogonal scores partial least squares 
#' regressions with the NIPALS algorithm. It return a comprehensive set of
#' pls outputs (e.g. scores and vip). 
#' @param X a O2pls object or a matrix of predictor variables.
#' @param Y a single vector indicate the group
#' @param nc the number of pls components (the one joint components + 
#'  number of orthogonal components ).
#' @param scale logical indicating whether \code{X} must be scaled (suggest TRUE).
#' @param center boolean values determining if data should be centered or not
#' @param maxiter maximum number of iterations.
#' @param tol limit for convergence of the algorithm in the nipals algorithm.
#' @return a list containing the following elements:
#' \itemize{
#' \item{\code{nc}}{ the number of components used(one joint components + 
#'  number of orthogonal components }
#' \item{\code{scores}}{ a matrix of scores corresponding to the observations 
#' in \code{X}, The components retrieved correspond to the ones optimized 
#' or specified.}
#' \item{\code{Xloadings}}{ a matrix of loadings corresponding to the
#'  explanatory variables. The components retrieved correspond to the ones
#'  optimized or specified.}
#' \item{\code{Yloadings}}{ a matrix of partial least squares loadings
#'  corresponding to \code{Y}}
#' \item{\code{vip}}{ the VIP matrix.}
#' \item{\code{xvar}}{ a matrix indicating the standard deviation of each
#'  component (sd), the variance explained by each single component
#'  (explained_var) and the cumulative explained variance
#'  (cumulative_explained_var). These values are
#'  computed based on the data used to create the projection matrices.}
#' \item{\code{projection_matrix}}{ the matrix of  projection matrix}
#' \item{\code{weight}}{ a matrix of partial least squares ("pls") weights.}} 
#' @examples 
#' X <- matrix(rnorm(50),10,5)
#' Y <- matrix(rnorm(50),10,5)
#' fit <- o2pls(X,Y,2,1,1)
#' yy <- rep(c(0,1),5)
#' fit0 <- oplsda(fit,yy,2)
#' @author Kai Guo
#' @export
oplsda <- function(X, Y, nc, scale=FALSE, center=TRUE, maxiter = 100, tol = 1E-5){
    if (inherits(x = X, what = "O2pls")) {
        result <- X@results
        X <- cbind(result$Xscore%*%t(result$Xloading),result$Yscore%*%t(result$Yloading))
    }
    X <- as.matrix(X)
    if(is.character(Y)){
        Y <- as.numeric(as.factor(Y))
    }
    Y <- as.matrix(cbind(Y))
    if(nrow(X)!=nrow(Y)) stop("Y should have same length as number of X rows")
    if(isTRUE(scale)){
        X = scale(X,center,scale=TRUE)
    }
    if(isTRUE(center)&!isTRUE(scale)){
        X = scale(X,center,scale=FALSE)
    }
    fit <- opls(X,Y,nc,maxiter,tol)
    score <- fit$scores
    colnames(score) <- paste0("LV",1:nc)
    rownames(score) <- rownames(X)
    Xloading <- t(fit$X_loadings)
    colnames(Xloading) <- paste0("LV",1:nc)
    rownames(Xloading) <- colnames(X)
    Yloading <- fit$Y_loadings
    vip <- fit$vip
    colnames(vip) <- paste0("LV",1:nc)
    rownames(vip) <- colnames(X)
    xvar <- fit$variance$x_var
    colnames(xvar) <- paste0("LV",1:nc)
    rownames(xvar) <- c("sd","explained_var","cumulative_explained_var")
    weight <- t(fit$weights)
    rownames(weight) <- colnames(X)
    colnames(weight) <- paste0("LV",1:nc)
    projection_matrix <- fit$projection_mat
    rownames(projection_matrix) <- colnames(X)
    colnames(projection_matrix) <- paste0("LV",1:nc)
    R2X <- s2(score%*%t(Xloading))/s2(X)
    res <- list(ncomp = nc, score = score, Xloading = Xloading,
              R2X = R2X,
              Yloading = Yloading,
              vip = vip, xvar = xvar,
              weight = weight,
              projection_matrix = projection_matrix)
    class(res) <- "o2plsda"
    return(res)
}