#' @title Partial least squares discriminant analysis
#' @description Perform a PLS discriminant analysis
#' @importFrom stats model.matrix cor 
#' @usage 
#' plsda(X, Y, nc, scale = TRUE, center = TRUE, cv = TRUE, nr_folds = 5)
#' @param X a matrix of predictor variables.
#' @param Y a single vector indicate the group
#' @param nc the number of pls components (the one joint components + 
#'  number of orthogonal components ).
#' @param scale logical indicating whether \code{X} must be scaled (suggest TRUE).
#' @param center logical indicating whether \code{X} must be centered (suggest TRUE).
#' @param cv logical indicating whether cross-validation will be performed or not (suggest TRUE).
#' @param nr_folds nr_folds Integer to indicate the folds for cross validation.
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
#' \item{\code{xvar}}{ variance explained of X by each single component.}
#' \item{\code{R2Y}}{ variance explained of Y by each single component.}
#' \item{code{PRESS}}{The residual sum of squares for the samples which were not used to fit the model}
#' \item{code{Q2}}{quality of  cross-validation}}
#' @examples 
#' X <- matrix(rnorm(500),10,50)
#' Y <- rep(c("a","b"),each=5)
#' fit <- plsda(X,Y,2)
#' @author Kai Guo
#' @export

plsda <- function(X,Y,nc,scale=TRUE, center=TRUE, cv = TRUE, nr_folds = 5){
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
    fit <- .pls(X,Y,nc,cv=cv,nr_folds=nr_folds)
    nc <- fit$nc
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
    R2Y <- fit$R2Y
    names(R2Y) <- paste0('LV',1:nc)
    res <- list(nc = nc, scores = score, Xloading = Xloading, vip = VIP, xvar = xvar, R2Y = R2Y, PRESS = fit$PRESS, Q2 = fit$Q2,Q2G = fit$Q2G)
    class(res)<-"plsda"
    return(res)
}

#' @title Partial least squares discriminant analysis
#' @description Perform a PLS discriminant analysis
#' @param X a matrix of predictor variables.
#' @param Y a single vector indicate the group
#' @param nc the number of pls components .
#' @param cv logical indicating whether cross-validation will be performed or not (suggest TRUE).
#' @param nr_folds nr_folds Integer to indicate the folds for cross validation.
#'  @return list with PLS results
#' @keywords internal
.pls <- function(X, Y , nc, cv = TRUE, nr_folds = 5){
    n = nrow(X)
    m = ncol(X)
    q = ncol(Y)
    px = ncol(X)
    Xt <- X
    Yt <- Y
    xvar <- sum((Xt)^2)
    x_var <- NULL
    if(isTRUE(cv)){
        Xs = svd(Xt, nu=0, nv=0)
        rank_X = sum(Xs$d > 0.0001)
        nc = min(n, rank_X)
    }
    if(nc == n) nc == n - 1
    W <- matrix(0, m, nc)
    U <- matrix(0, n, nc)
    Tt <- matrix(0, n, nc)
    Ch <- matrix(0, q, nc)
    Ph <- matrix(0, px, nc)
    bh <- rep(0, nc)
    RSS = rbind(rep(n-1,q), matrix(NA, nc, q))
    PRESS = matrix(NA, nc, q)
    Q2 = matrix(NA, nc, q)
    fold = split(sample(1:n), rep(1:nr_folds, length=n))
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
        #### cross validation modified from 
        #### https://github.com/gastonstat/DiscriMiner/blob/master/R/my_plsDA.R
        RSS[i+1,] = colSums((Yt - tp%*%t(ct))^2)
        press = matrix(0, n, q)
        if(isTRUE(cv)){
           for (k in 1:nr_folds){
               omit=fold[[k]]
               uh.cv <- Yt[-omit,1]
               wh.cvt <- rep(1,px)
                itcv <- 1
                repeat{
                    wh.cv <- t(Xt[-omit,]) %*% uh.cv / sum(uh.cv^2)
                    wh.cv <- wh.cv / sqrt(sum(wh.cv^2))
                    th.cv <- Xt[-omit,] %*% wh.cv
                    ch.cv <- t(Yt[-omit,]) %*% th.cv / sum(th.cv^2)
                    uh.cv <- Yt[-omit,] %*% ch.cv / sum(ch.cv^2)
                    wcv.dif <- wh.cv - wh.cvt
                    wh.cvt <- wh.cv
                    if (sum(wcv.dif^2)<1e-06 || itcv==100) break
                    itcv <- itcv + 1
                    }
                    Yhat.cv = (Xt[omit,] %*% wh.cv) %*% t(ch.cv)
                    press[omit,] = (Yt[omit,] - Yhat.cv)^2
           }
                    PRESS[i,] = colSums(press)
                    Q2[i,] = 1 - PRESS[i,]/RSS[i,]
                    }
        # deflation
        Xt = Xt - (tp %*% t(p))
        Yt = Yt - (tp %*% t(ct))
        W[,i] = w
        U[,i] = u
        Tt[,i] = tp
        Ch[,i] = ct
        Ph[,i] = p
        bh[i] = t(u) %*% tp
        Q2G = 1 - rowSums(PRESS)/rowSums(RSS[-nc,])
    }
    ## modified from the mixOmics package
    if(isTRUE(cv)){
        ns = max(which(Q2G >= 0.05))
        if(ns < nc){
            nc <- ns
            if(nc < 2){
                nc <- 2
            }
            W = W[,1:nc, drop=FALSE]
            U = U[,1:nc, drop=FALSE]
            Ph = Ph[,1:nc,drop=FALSE]
            Tt = Tt[,1:nc,drop=FALSE]
            Ch = Ch[,1:nc,drop=FALSE]
        }
    }
    for (j in 1:nc)
    {
        v <- t(Tt[, j, drop=FALSE]) %*% X
        varx <- v%*%t(v) /crossprod(Tt[, j],Tt[, j])/xvar
        x_var = c(x_var, varx)
    }
    ### https://github.com/gastonstat/DiscriMiner/blob/master/R/my_plsDA.R
    ty = cor(Y, Tt)
    R2y = cor(Y, Tt)^2 # R2 coefficients
    R2Y = colMeans(R2y) # Sum of squares Y
    return(list(nc = nc, W = W,U = U,Tt = Tt,Ch = Ch, Ph = Ph, x_var = x_var, R2Y = R2Y, PRESS = PRESS, Q2 = Q2, Q2G = Q2G))
}



