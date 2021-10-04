# OPLSDA - where X is the concatenated joint variance projection: cbind(Tp %*% t(Wp), Up %*% t(Cp)), 
# y is the variable to regress and n is the desired number of components
#' @ still under developing
#' @export
oplsda <- function(fit, y, n) {
    if(is.character(y)){
        y<-as.numeric(as.factor(y))
    }
    X <- cbind(fit$Tp %*% t(fit$Wp), fit$Up %*% t(fit$Cp))
    To <- Wo <- Po <- NULL
    k <- 1
    while(k <= n) {
        w  <- t(t(y) %*% X / sum(y^2))           # 1
        w  <- w / sqrt(sum(w^2))                 # 2
        tt <- X %*% w                            # 3
        cc <- t( t(tt) %*% y / sum(tt^2)  )      # 4
        u  <- y %*% cc / sum(cc^2)               # 5
        p  <- t( t(tt) %*% X / sum(tt^2)  )      # 6
        wo <- p - sum(w * p) * w                 # 7
        wo <- wo / sqrt(sum(wo^2))               # 8
        tto <- X %*% wo                          # 9
        po  <- t( t(tto) %*% X / sum(tto^2) )    # 10
        Eopls  <- X - tto %*% t(po)              # 11
        To <- cbind(To, tto)
        Po <- cbind(Po, po)
        Wo <- cbind(Wo, wo)
        
        X <- Eopls 
        k <- k+1
    }
    list(X=X, To=To, Po=Po, Wo=Wo)
}