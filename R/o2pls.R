#' fit O2PLS model with best nc, nx, ny
#' @importFrom methods new
#' @param X a Numeric matrix (input)
#' @param Y a Numeric matrix (input)
#' @param nc Integer. Number of joint PLS components.
#' @param nx Integer. Number of orthogonal components in X 
#' @param ny Integer. Number of orthogonal components in Y 
#' @param scale boolean values determining if data should be scaled or not 
#' @param center boolean values determining if data should be centered or not
#' @return An object containing
#'    \item{Xscore}{Joint \eqn{X} scores}
#'    \item{Xloading}{Joint \eqn{X} loadings}
#'    \item{Yscore}{Joint \eqn{Y} scores}
#'    \item{Yloading}{Joint \eqn{Y} loadings}
#'    \item{TYosc}{Orthogonal \eqn{X} scores}
#'    \item{PYosc}{Orthogonal \eqn{X} loadings}
#'    \item{WYosc}{Orthogonal \eqn{X} weights}
#'    \item{UXosc}{Orthogonal \eqn{Y} scores}
#'    \item{PXosc}{Orthogonal \eqn{Y} loadings}
#'    \item{CXosc}{Orthogonal \eqn{Y} weights}
#'    \item{BU}{Regression coefficient in \code{Tt} ~ \code{U}}
#'    \item{BT}{Regression coefficient in \code{U} ~ \code{Tt}}
#'    \item{R2Xhat}{Prediction of \eqn{X} with \eqn{Y}}
#'    \item{R2Yhat}{Prediction of \eqn{Y} with \eqn{X}}
#'    \item{R2X}{Variation of the modeled part in \eqn{X} (defined by Joint + Orthogonal variation) as proportion of total variation in \eqn{X}}
#'    \item{R2Y}{Variation of the modeled part in \eqn{Y} (defined by Joint + Orthogonal variation) as proportion of total variation in \eqn{Y}}
#'    \item{R2Xcorr}{Variation of the joint part in \eqn{X} }
#'    \item{R2Ycorr}{Variation of the joint part in \eqn{Y} }
#'    \item{R2Xo}{Variation of the orthogonal part in \eqn{X} as proportion of variation in \eqn{X}}
#'    \item{R2Yo}{Variation of the orthogonal part in \eqn{Y} as proportion of variation in \eqn{Y}}
#'    \item{R2Xp}{Variation in \eqn{X} joint part predicted by \eqn{Y} Joint part}
#'    \item{R2Yp}{Variation in \eqn{Y} joint part predicted by \eqn{X} Joint part}
#'    \item{varXj}{Variation in each Latent Variable (LV) in \eqn{X} Joint part}
#'    \item{varYj}{Variation in each Latent Variable (LV) in \eqn{Y} Joint part}
#'    \item{varXorth}{Variation in each Latent Variable (LV) in \eqn{X} Orthogonal part}
#'    \item{varYorth}{Variation in each Latent Variable (LV) in \eqn{Y} Orthogonal part}
#'    \item{Exy}{Residuals in \eqn{X}}
#'    \item{Fxy}{Residuals in \eqn{Y}}
#'    
#' @examples 
#' set.seed(123)
#' X = matrix(rnorm(500),50,10)
#' Y = matrix(rnorm(500),50,10)
#' X = scale(X, scale = TRUE)
#' Y = scale(Y, scale = TRUE)
#' fit <- o2pls(X, Y, 1, 2, 2)
#' summary(fit)
#' @export
#' @author Kai Guo

o2pls <- function(X,Y,nc,nx,ny,scale=FALSE,center=FALSE){
    ## total variances
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    if(isTRUE(scale)){
        X = scale(X,center,scale=TRUE)
        Y = scale(Y,center,scale=TRUE)
    }
    if(isTRUE(center)&!isTRUE(scale)){
        X = scale(X,center,scale=FALSE)
        Y = scale(Y,center,scale=FALSE)
    }
    if(ncol(X) < max(nc)+max(nx,ny) | ncol(Y) < max(nc)+max(nx,ny))
        message("The combinations of # components should be fewer than the data dimensions\n")
    if(ncol(X) < min(nc)+min(nx,ny) | ncol(Y) < min(nc)+min(ny,ny))
        stop("There is no valid combination of numbers of components! Please select fewer components in nc, nx, ny.\n")
    if (nc != round(abs(nc)) || nc <= 0) {
        stop("nc should be a positive integer\n")
    }
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
    Xt <- X
    Yt <- Y
    n <- nrow(X)
    ncx<- ncol(X)
    ncy <- ncol(Y)
    SSX <- s2(X)
    SSY <- s2(Y)
    #define orth
    TYosc <- matrix(0, n, nnx)
    PYosc <- matrix(0, ncx, nnx)
    WYosc <- matrix(0, ncx, nnx)
    UXosc <- matrix(0, n, nny)
    PXosc <- matrix(0, ncy, nny)
    CXosc <- matrix(0, ncy, nny)
    
    # Calculate the PCA components of Y'T
   ## nn<-nc+max(nx+ny)
    CDW <- eigenmult(t(Y),X)
    cdw <- svd(CDW,nu=nc,nv=nc)
    Xloading <- cdw$v
    Yloading <- cdw$u
   # Xscore <- X %*% Xloading
    Xscore <- eigenmult(X,Xloading)
#    Yscore <- Y %*% Yloading
    Yscore <- eigenmult(Y,Yloading)
    if(nx>0){
        for(i in 1:nx){
          #  Exy <- X - Xscore%*%t(Xloading)
            Exy <- X - eigenmult(Xscore, t(Xloading))
          #  wosc<- svd(t(Exy)%*%Xscore,nu=1,nv=0)$u
            wosc<- svd(eigenmult(t(Exy),Xscore),nu=1,nv=0)$u
           # t_yo <- X%*%wosc
            t_yo <- eigenmult(X,wosc)
         #   p_yo <- t(solve(t(t_yo)%*%t_yo)%*%t(t_yo)%*%X)
            p_yo <- t(eigenthree(solve(AtA(t_yo)),t(t_yo),X))
          #  X <- X - t_yo%*%t(p_yo)
            X <- X - eigenmult(t_yo,t(p_yo))
            TYosc[,i] <- t_yo
            PYosc[,i] <- p_yo
            WYosc[,i] <- wosc	
            #Xscore <- X%*%Xloading
            Xscore <- eigenmult(X, Xloading)
        }
    }
    if(ny>0){
        for(i in 1:ny){
           # Fxy <- Y - Yscore %*% t(Yloading)
            Fxy <- Y - eigenmult(Yscore, t(Yloading))
           # cosc <- svd(t(Fxy)%*%Yscore,nu=1,nv=0)$u
            cosc <- svd(eigenmult(t(Fxy), Yscore),nu=1,nv=0)$u
           # u_xo <- Y%*%cosc
            u_xo <- eigenmult(Y, cosc)
           # p_xo <- t(solve(t(u_xo)%*%u_xo)%*%t(u_xo)%*%Y)
            p_xo <- t(eigenthree(solve(AtA(u_xo)),t(u_xo),Y))
          #  Y <- Y - u_xo%*%t(p_xo)
            Y <- Y - eigenmult(u_xo, t(p_xo))
            UXosc[,i] <- u_xo
            PXosc[,i] <- p_xo
            CXosc[,i] <- cosc	
          #  Yscore <- Y%*%Yloading
            Yscore <- eigenmult(Y, Yloading)
        }
    }
    BU <- eigenthree(solve(AtA(Yscore)),t(Yscore),Xscore)
    BT <- eigenthree(solve(AtA(Xscore)),t(Xscore),Yscore)
    Exy <- Xt - eigenmult(Xscore, t(Xloading)) - eigenmult(TYosc, t(PYosc))
    Fxy <- Yt - eigenmult(Yscore, t(Yloading)) - eigenmult(UXosc, t(PXosc))
    R2Xcorr <- s2(eigenmult(Xscore, t(Xloading)))/SSX
    R2Ycorr <- s2(eigenmult(Yscore, t(Yloading)))/SSY
    R2Xo <- s2(eigenmult(TYosc, t(PYosc)))/SSX
    R2Yo <- s2(eigenmult(UXosc, t(PXosc)))/SSY
    R2X <- 1- s2(Exy)/SSX
    R2Y <- 1- s2(Fxy)/SSY
    Xh <- eigenthree(Yscore, BU, t(Xloading))
    Yh <- eigenthree(Xscore, BT, t(Yloading))
    R2Xhat <- 1 - s2(Xh-Xt)/SSX
    R2Yhat <- 1 - s2(Yh-Yt)/SSY
    R2Xp <- (s2(Xh)/SSX)/R2Xcorr
    R2Yp <- (s2(Yh)/SSY)/R2Ycorr
    Qx <- Q(Xt,Xh)
    Qy <- Q(Yt,Yh)
    varXj = apply(Xscore,2,function(x)sum(x^2))/SSX
    varYj = apply(Yscore,2,function(x)sum(x^2))/SSY
    varXorth = apply(PYosc,2,function(x)sum(x^2))*apply(TYosc,2,function(x)sum(x^2))/SSX
    varYorth = apply(PXosc,2,function(x)sum(x^2))*apply(UXosc,2,function(x)sum(x^2))/SSY
    rownames(Xscore) <- rownames(TYosc) <- rownames(Exy) <- rownames(Xt)
    rownames(Yscore) <- rownames(UXosc) <- rownames(Fxy) <- rownames(Yt)
    rownames(Xloading) <- rownames(PYosc) <- rownames(WYosc) <- colnames(Exy) <- colnames(Xt)
    rownames(Yloading) <- rownames(PXosc) <- rownames(CXosc) <- colnames(Fxy) <- colnames(Yt)
    colnames(Xscore) <- colnames(Yscore) <- colnames(Xloading) <- colnames(Yloading) <- paste0("LV",1:nc)
    colnames(TYosc) <- colnames(PYosc) <- colnames(WYosc) <- paste0("LV", 1:nnx)
    colnames(UXosc) <- colnames(PXosc) <- colnames(CXosc) <- paste0("LV", 1:nny)
    res <- list(Xscore=Xscore,Yscore=Yscore,Xloading=Xloading,Yloading=Yloading,
                R2Xcorr=R2Xcorr,R2Ycorr=R2Ycorr,
                BU = BU, BT = BT,
                R2Xo = R2Xo, R2Yo = R2Yo,
                R2Xp = R2Xp, R2Yp = R2Yp,
                TYosc = TYosc, UXosc = UXosc,
                PXosc = PXosc, PYosc = PYosc,
                WYosc = WYosc, CXosc = CXosc,
                R2Xhat = R2Xhat, R2Yhat = R2Yhat,
                Qx = Qx, Qy = Qy,
                varXj = varXj, varYj = varYj, 
                varXorth = varXorth, varYorth =varYorth,
                Exy = Exy, Fxy =Fxy,
                R2X = R2X, R2Y = R2Y)
    results <- new("O2pls",
                   X = Xt,
                   Y = Yt,
                   params = list(nc = nc, nx = nx, ny = ny),
                   results = res)
    return(results)
}
