#' @title Cross validation for O2PLS
#' @importFrom dplyr group_by arrange desc
#' @importFrom dplyr summarise mutate filter
#' @importFrom magrittr %>%
#' @importFrom parallel makePSOCKcluster clusterEvalQ clusterExport parLapply
#' @importFrom parallel mclapply stopCluster
#' @param X a Numeric matrix (input)
#' @param Y a Numeric matrix (input)
#' @param nc Integer. Number of joint PLS components.
#' @param nx Integer. Number of orthogonal components in X 
#' @param ny Integer. Number of orthogonal components in Y 
#' @param group a vector to indicate the group for Y
#' @param nr_folds Integer to indicate the folds for cross validation
#' @param ncores Integer. Number of CPUs to use for cross validation
#' @param scale boolean values determining if data should be scaled or not 
#' @param center boolean values determining if data should be centered or not
#' @examples 
#' set.seed(123)
#' X = matrix(rnorm(500),50,10)
#' Y = matrix(rnorm(500),50,10)
#' X = scale(X, scale = TRUE)
#' Y = scale(Y, scale = TRUE)
#' # group factor could be omitted if you don't have any group 
#' group <- rep(c("Ctrl","Treat"), each = 25)
#' cv <- o2cv(X, Y, 1:2, 1:2, 1:2, group=group, nr_folds = 2, ncores=1)
#' @author Kai Guo
#' @return a data frame with the Q and RMSE values
#' @export
#' 
o2cv<-function(X, Y, nc, nx, ny, group=NULL, nr_folds = 5, ncores=1,
            scale = FALSE, center = FALSE){
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
    if(nrow(X) < nr_folds) stop("There are more folds than samples, please set nr_folds <= ",nrow(X),"\n")
    stopifnot(ncores == abs(round(ncores)))
    if(nr_folds == 1){stop("At least two folds needed in the Cross-validation \n")}
    ###build params data.frame
    if(is.null(group)){
        group <- rep(1,nrow(X))
    }else{
        if(length(group)!=nrow(X)) stop("The group should have same length of sample\n")
    }
    params <- data.frame(t(expand.grid(nc,nx,ny)))
    colnames(params)<-NULL
    rownames(params)<-c("n","nx","ny")
    cl <- NULL
    on.exit({if(!is.null(cl)) stopCluster(cl)})
    if(Sys.info()[["sysname"]] == "Windows" && ncores > 1){
        cl <- makePSOCKcluster(ncores)
        clusterEvalQ(cl, library(o2plsda))
        clusterExport(cl, varlist = ls(), envir = environment())
        res <- parLapply(cl,params,function(f){
            suppressMessages(.o2cv(X,Y,n=f[1],nx=f[2],ny=f[3],group,nr_folds))
        })
    } else {
        res <- mclapply(mc.cores = ncores,params,function(f){
            .o2cv(X,Y,n=f[1],nx=f[2],ny=f[3],group,nr_folds)
        })
    }
    res <- do.call("c",res)
    results <- as.data.frame(t(sapply(res, function(x)unlist(x))))
    results <- results%>%group_by(nc,nx,ny)%>%
        summarise(Qx=mean(Qx),Qy=mean(Qy),Px=mean(Px),Py=mean(Py),Rx=mean(Rx),Ry=mean(Ry))%>%
        mutate(RMSE=Rx+Ry,Qxy=sqrt(Qx^2+Qy^2))%>%arrange(desc(Qxy))%>%filter(Qx>=0,Qy>=0)
    nc <- as.data.frame(results)[1,1]
    nx <- as.data.frame(results)[1,2]
    ny <- as.data.frame(results)[1,3]
    Qxy <- as.data.frame(results)[1,11]
    RMSE <- as.data.frame(results)[1,10]
    message("#####################################")
    message("The best parameters are nc = ",nc,", nx = ",nx,", ny = ",ny)
    message("#####################################")
    message("The Qxy is ",Qxy, " and the RMSE is: ", RMSE)
    message("#####################################")
    return(results)
}


#' @title do cross-validation with group factors
#' @param X a Numeric matrix (input)
#' @param Y a Numeric matrix (input)
#' @param n Integer. Number of joint PLS components.
#' @param nx Integer. Number of orthogonal components in X 
#' @param ny Integer. Number of orthogonal components in Y 
#' @param group a vector to indicate the group for Y
#' @param nr Integer to indicate the folds for cross validation
#' @return a list of o2pls results
#' @keywords internal
#' @author Kai Guo

.o2cv<-function(X,Y,n,nx,ny,group,nr){
    results <- list()
    cls.grp <- getMCCV_cpp(group, n = nr)
    for(k in 1:max(cls.grp,na.rm=TRUE)) {
        # selects training group
        Xtr <- X[cls.grp != k & is.na(cls.grp)==FALSE,]
        Ytr <- Y[cls.grp != k & is.na(cls.grp)==FALSE,]
        #do o2pls with training data
        o2 <- o2pls(Xtr, Ytr, n, nx ,ny)
        o2 <- o2@results
        #subset testing data
        Xev <- X[cls.grp == k & is.na(cls.grp)==FALSE,]
       
        
        # Predicts y_hat
        if(nx>0){
            To <- matrix(0,nrow(Xev),nx)
        for (i in 1:nx) {
            #tx <- Xev %*% o2$WYosc[ , i,drop=F]
            tx <- eigenmult(Xev, o2$WYosc[ , i,drop=F])
          #  Xev <- Xev - tx %*% t(o2$PYosc[ , i,drop=F])
            Xev <- Xev - eigenmult(tx, t(o2$PYosc[ , i,drop=F]))
            To[,i]<- tx
        }}
      #  Tpp <- Xev %*% o2$Xloading
      #  Y_hat <- Tpp %*% o2$BT %*% t(o2$Yloading)
        Tpp <- eigenmult(Xev, o2$Xloading)
        Y_hat <- eigenthree(Tpp,o2$BT, t(o2$Yloading))
        
        # Predicts x_hat
        Yev <- Y[cls.grp == k & is.na(cls.grp)==FALSE,]
        
        if(ny>0){
            Uo <- matrix(0,nrow(Yev),ny)
        for (i in 1:ny) {
          #  ux <-  Yev %*% o2$CXosc[ , i,drop=F]
            ux <-  eigenmult(Yev, o2$CXosc[ , i,drop=F])
           # Yev <- Yev - ux %*% t(o2$PXosc[ , i,drop=F])
            Yev <- Yev - eigenmult(ux, t(o2$PXosc[ , i,drop=F]))
            Uo[,i] <- ux
        }
        }
     #   Upp <- Yev %*% o2$Yloading
     #   X_hat <- Upp %*% o2$BU %*% t(o2$Xloading)
        Upp <- eigenmult(Yev, o2$Yloading)
        X_hat <- eigenthree(Upp, o2$BU, t(o2$Xloading))
        
        # Restore to original values
        Xev <- X[cls.grp == k & is.na(cls.grp)==FALSE,]
        Yev <- Y[cls.grp == k & is.na(cls.grp)==FALSE,]
        
        tmp <- list(k=k, nc=n, nx=nx, ny=ny, Qx=Q(Xev,X_hat),
                    Qy=Q(Yev,Y_hat),Px=s2(Xev-X_hat),Py=s2(Yev-Y_hat),
                    Rx=rcpp_rmse(Xev,X_hat),Ry=rcpp_rmse(Yev,Y_hat))
        results <- append(results, list(tmp))
    }
    return(results)
}
