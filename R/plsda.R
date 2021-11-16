#' @title Partial least squares discriminant analysis
#' @description Perform a PLS discriminant analysis
#' @importFrom stats model.matrix cor 
#' @usage 
#' plsda(X, Y, nc, scale, center)
#' @param X a matrix of predictor variables.
#' @param Y a single vector indicate the group
#' @param nc the number of pls components (the one joint components + 
#'  number of orthogonal components ).
#' @param scale logical indicating whether \code{X} must be scaled (suggest TRUE).
#' @param center logical indicating whether \code{X} must be centered (suggest TRUE).
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
#' \item{\code{vip}}{ the VIP matrix.}
#' \item{\code{xvar}}{ variance explained by each single component}}
#' @examples 
#' X <- matrix(rnorm(50),10,5)
#' Y <- rep(c(0,1),each=5)
#' fit <- plsda(X,Y,2)
#' @author Kai Guo
#' @export

plsda <- function(X,Y,nc,scale=TRUE, center=TRUE){
    X <- as.matrix(X)
    Y <- model.matrix(~ -1 + Y)
    if(nrow(X)!=nrow(Y)) stop("Y should have same length as number of X rows")
    if(isTRUE(scale)){
        X = scale(X,center,scale=TRUE)
        Y = scale(Y,center,scale=TRUE)
    }
    if(isTRUE(center)&!isTRUE(scale)){
        X = scale(X,center,scale=FALSE)
        Y = scale(Y,center,scale=FALSE)
    }
    fit <- .pls(X,Y,nc)
    score <- fit$Tt
    colnames(score) <- paste0("LV",1:nc)
    rownames(score) <- rownames(X)
    Xloading <- fit$W
    colnames(Xloading) <- paste0("LV",1:nc)
    rownames(Xloading) <- colnames(X)
    ## calculate the VIP values
    W <- Xloading
    H <- nc
    q <- ncol(Y)
    p <- ncol(X)
    Th <- score
    VIP <- matrix(0, nrow = p, ncol = H)
    cor2 <- cor(Y, Th, use = "pairwise")^2
    cor2 <- as.matrix(cor2, nrow = q)
    VIP[, 1] <- W[, 1]^2
    if (H > 1) {
        for (h in 2:H) {
            if (q == 1) {
                R = cor2[, 1:h]
                VIP[, h] = R %*% t(W[, 1:h]^2)/sum(R)
            }
            else {
                R = apply(cor2[, 1:h], 2, sum)
                VIP[, h] = R %*% t(W[, 1:h]^2)/sum(R)
            }
        }
    }
    VIP <- sqrt(p * VIP)
    rownames(VIP) <- rownames(W)
    colnames(VIP) <- paste0("comp", 1:H)
    xvar <- fit$x_var
    names(xvar) <- paste0('LV',1:nc)
    res <- list(nc = nc, scores = score, Xloading = Xloading, vip = VIP, xvar = xvar)
    class(res)<-"plsda"
    return(res)
}

#' @title Partial least squares discriminant analysis
#' @description Perform a PLS discriminant analysis
#' @param X a matrix of predictor variables.
#' @param Y a single vector indicate the group
#' @param nc the number of pls components (the one joint components + 
#'  number of orthogonal components ).
#'  @return list with PLS results
#' @keywords internal
.pls <- function(X,Y,nc){
    n = nrow(X)
    m = ncol(X)
    q = ncol(Y)
    Xt <- X
    Yt <- Y
    xvar <- sum((Xt)^2)
    x_var <- NULL
    W <- matrix(0, m, nc)
    U <- matrix(0, n, nc)
    Tt <- matrix(0, n, nc)
    Ch <- matrix(0, q, nc)
    bh <- rep(0, nc)
    ## PLS2 algorithm
    for (i in 1:nc)
    {
        u = Yt[,1]
        wo = rep(1, m)
        iter = 1
        repeat
        {
            w = t(Xt) %*% u / sum(u^2)
            w = w / sqrt(sum(w^2))# normalize wo
            tp = Xt %*% w
            ct = t(Yt) %*% tp / sum(tp^2)
            u = Yt %*% ct / sum(ct^2)
            w.dif = w - wo
            wo = w
            if (sum(w.dif^2)<1e-06 || iter==100) break
            iter = iter + 1
        }
        p = t(Xt) %*% tp / sum(tp^2)
        # deflation
        Xt = Xt - (tp %*% t(p))
        Yt = Yt - (tp %*% t(ct))
        W[,i] = w
        U[,i] = u
        Tt[,i] = tp
        Ch[,i] = ct
        bh[i] = t(u) %*% tp
    }
    ## modified from the mixOmics package
    for (j in 1:nc)
    {
        v <- t(Tt[, j, drop=FALSE]) %*% X
        varx <- v%*%t(v) /crossprod(Tt[, j],Tt[, j])/xvar
        x_var = c(x_var, varx)
    }
    return(list(W=W,U=U,Tt=Tt,Ch=Ch,x_var=x_var))
}



