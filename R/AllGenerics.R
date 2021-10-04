#' Extract the loadings from an O2PLS fit
#'
#' This function extracts loading parameters from an O2PLS fit
#'
#' @param x Object of class \code{O2pls}
#' @param ... For consistency
#' 
#' @return Loading matrix
#' 
#' @rdname loadings
#' @export
loadings <- function(x, ...) UseMethod("loadings")

#' @export
loadings.O2pls <- function(x,loading = c("Xjoint", "Yjoint", "Xorth", "Yorth"),...){
    x<-x@results
    if(loading=="Xjoint"){
        res <- x$Xloading
    }else if(loading == "Yjoint"){
        res <- x$Yloading
    }else if(loading == "Xorth"){
        res <- x$PYosc
    }else if(loading == "Yorth"){
        res <- x$PXosc
    }else{
        stop('Please specify the loading: ["Xjoint", "Yjoint", "Xorth", "Yorth"] \n')
    }
    return(res)
    
}
#' Extract the scores from an O2PLS fit
#'
#' This function extracts score matrices from an O2PLS fit
#'
#' @param x Object of class \code{O2pls}
#' @param ... For consistency
#' 
#' @return Scores matrix
#' 
#' 
#' @rdname scores
#' @export
scores <- function(x, ...) UseMethod("scores")

#' @export
scores.O2pls <- function(x, score = c("Xjoint", "Yjoint", "Xorth", "Yorth"),...){
    x<-x@results
    if(score=="Xjoint"){
        res <- x$Xscore
    }else if(score == "Yjoint"){
        res <- x$Yscore
    }else if(score == "Xorth"){
        res <- x$TYosc
    }else if(score == "Yorth"){
        res <- x$UXosc
    }else{
        stop('Please specify the score: ["Xjoint", "Yjoint", "Xorth", "Yorth"] \n')
    }
    return(res)
    
}