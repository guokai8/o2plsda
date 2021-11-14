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
    colnames(varj) <- paste0("LV",1:object@params$nc)
    rownames(varj) <- c("X","Y")
    varx <- rbind(res$varXorth)
    colnames(varx) <- paste0("LV",1:object@params$nx)
    rownames(varx) <- "X"
    vary <- rbind(res$varYorth)
    colnames(vary) <- paste0("LV",1:object@params$ny)
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

#' Extract the VIP values from the O2PLS-DA object
#' @param x the o2plsda object or plsda object
#' @return a data frame 
#' @export
vip <- function(x){
    res <- x$vip
    return(res)
    
}
