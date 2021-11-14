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
#' Extract the loadings from an O2PLS fit
#'
#' This function extracts loading parameters from an O2PLS fit
#'
#' @param x Object of class \code{O2pls}
#' @param loading the loadings for one of "Xjoint", "Yjoint", "Xorth", "Yorth"
#' @param ... For consistency
#' @return Loading matrix
#' 
#' @rdname loadings
#' @export
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
#' @title extract the loading value from the O2PLSDA analysis
#' @param x Object of class \code{o2plsda}
#' @param loading the loadings for one of "Xjoint", "Yjoint", "Xorth", "Yorth"
#' @param ... For consistency
#' @export
loadings.o2plsda <- function(x,loading ="Xloading",...){
    if(loading=="Xloading"){
        res <- x$Xloading
    }else if(loading == "Yloading"){
        res <- x$Yloading
    }else{
        stop('Please specify the loading: ["Xjoint", "Yjoint", "Xorth", "Yorth"] \n')
    }
    return(res)
    
}

#' @title extract the loading value from the PLSDA analysis
#' @param x Object of class \code{plsda}
#' @param ... For consistency
#' @export
loadings.plsda <- function(x,...){
    res <- x$Xloading
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

#' Extract the scores from an O2PLS fit
#'
#' This function extracts scores parameters from an O2PLS fit
#'
#' @param x Object of class \code{O2pls}
#' @param score the scores matrix for one of "Xjoint", "Yjoint", "Xorth", "Yorth"
#' @param ... Other arguments 
#' @return score matrix
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
#' Extract the scores from an O2PLS DA analysis
#'
#'
#' @param x Object of class \code{o2plsda}
#' @param ... Other arguments 
#' @return score matrix
#' @export
#' @author Kai Guo
scores.o2plsda <- function(x,...){
    res <- x$score
    return(res)
    
}
#' Extract the scores PLSDA analysis
#'
#'
#' @param x Object of class \code{plsda}
#' @param ... Other arguments 
#' @return score matrix
#' @export
#' @author Kai Guo
scores.plsda <- function(x,...){
    res <- x$score
    return(res)
    
}
