#' @title Summary of an O2PLS fit
#'
#' @param fit a O2pls object
#' @return Detail of O2PLS results
#' @author Kai Guo
#' @export
summary.O2pls<-function(fit){
    cat("\n######### Summary of the O2PLS results #########\n")
    cat("### Call o2pls(X, Y, nc=",fit@params$nc,", nx=",fit@params$nx,", ny=",fit@params$ny,") ###\n")
    sx <- s2(fit@X)
    sy <- s2(fit@Y)
    res <- fit@results
    d <-data.frame(X=c(res$R2Xcorr,res$R2Xo,1-res$R2X),Y=c(res$R2Ycorr,res$R2Yo,1-res$R2Y))
    varj <- as.data.frame(rbind(res$varXj,res$varYj))
    colnames(varj) <- paste0("LV",1:fit@params$nc)
    rownames(varj) <- c("X","Y")
    varx <- rbind(res$varXorth)
    colnames(varx) <- paste0("LV",1:fit@params$nx)
    rownames(varx) <- "X"
    vary <- rbind(res$varYorth)
    colnames(vary) <- paste0("LV",1:fit@params$ny)
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
#' @return NULL
#' @author Kai Guo
#' @export
print.O2pls <- function (fit, ...) {
    summary(fit)
}


#' @export
vip <- function(x,...){
    res <- x$vip
    return(res)
    
}
