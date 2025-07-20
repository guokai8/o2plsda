
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
#' Extract selected variables from sparse models
#' @export
selected_vars <- function(x, ...) {
    UseMethod("selected_vars")
}
#' Extract selected variables from sparse PLS-DA models
#' 
#' @param x A sparse_plsda object
#' @param ... Additional arguments (currently unused)
#' @return Vector of selected variable indices
#' @export
#' @export
selected_vars.sparse_plsda <- function(x, ...) {
    if(inherits(x, "sparse_plsda")) {
        return(x@selected_vars)
    } else {
        return(x$selected_vars)
    }
}

#' Extract selected variables from SparseO2pls S4 objects
#' 
#' @param x A SparseO2pls S4 object
#' @param type Character. Which variables to extract: "X", "Y", or "both"
#' @param ... Additional arguments (currently unused)
#' @return Vector of variable indices (for "X" or "Y") or list of indices (for "both")
#' @export
selected_vars.SparseO2pls <- function(x, type = "X", ...) {
    if(!inherits(x, "SparseO2pls")) {
        stop("Object must be of class 'SparseO2pls'")
    }
    
    if(type == "X") {
        return(x@sparsity$selected_X)
    } else if(type == "Y") {
        return(x@sparsity$selected_Y)
    } else if(type == "both") {
        return(list(X = x@sparsity$selected_X, Y = x@sparsity$selected_Y))
    } else {
        stop("type must be 'X', 'Y', or 'both'")
    }
}

#' Default method for selected_vars
#' 
#' @param x An object
#' @param ... Additional arguments
#' @export
selected_vars.default <- function(x, ...) {
    stop("selected_vars method not implemented for objects of class: ", 
         paste(class(x), collapse = ", "))
}

#' @export
selected_vars.sparse_o2pls <- function(x, type = "X", ...) {
    # Handle both S4 and list objects
    if(inherits(x, "SparseO2pls")) {
        sparsity_info <- x@sparsity
    } else if(!is.null(x$sparsity)) {
        sparsity_info <- x$sparsity
    } else {
        stop("Cannot find sparsity information in object")
    }
    
    if(type == "X") {
        return(sparsity_info$selected_X)
    } else if(type == "Y") {
        return(sparsity_info$selected_Y)
    } else if(type == "both") {
        return(list(X = sparsity_info$selected_X, Y = sparsity_info$selected_Y))
    } else {
        stop("type must be 'X', 'Y', or 'both'")
    }
}

#' Get names of selected variables
#' 
#' @param x A sparse model object
#' @param type Character. Which variables to extract: "X", "Y", or "both"
#' @param ... Additional arguments passed to selected_vars
#' @return Character vector of variable names or list of name vectors
#' @export
selected_var_names <- function(x, type = "X", ...) {
    UseMethod("selected_var_names")
}

#' @export
selected_var_names.sparse_o2pls <- function(x, type = "X", ...) {
    indices <- selected_vars(x, type = type, ...)
    
    if(type == "both") {
        X_names <- if(!is.null(colnames(x$X))) colnames(x$X)[indices$X] else paste0("X", indices$X)
        Y_names <- if(!is.null(colnames(x$Y))) colnames(x$Y)[indices$Y] else paste0("Y", indices$Y)
        return(list(X = X_names, Y = Y_names))
    } else if(type == "X") {
        return(if(!is.null(colnames(x$X))) colnames(x$X)[indices] else paste0("X", indices))
    } else if(type == "Y") {
        return(if(!is.null(colnames(x$Y))) colnames(x$Y)[indices] else paste0("Y", indices))
    }
}

#' @export
selected_var_names.sparse_plsda <- function(x, ...) {
    indices <- selected_vars(x, ...)
    if(!is.null(x$variable_names)) {
        return(x$variable_names[indices])
    } else if(!is.null(rownames(x$loadings))) {
        return(rownames(x$loadings)[indices])
    } else {
        return(paste0("Var", indices))
    }
}

#' @export
selected_var_names.SparseO2pls <- function(x, type = "X", ...) {
    indices <- selected_vars(x, type = type, ...)
    
    if(type == "both") {
        X_names <- if(!is.null(colnames(x@X))) colnames(x@X)[indices$X] else paste0("X", indices$X)
        Y_names <- if(!is.null(colnames(x@Y))) colnames(x@Y)[indices$Y] else paste0("Y", indices$Y)
        return(list(X = X_names, Y = Y_names))
    } else if(type == "X") {
        return(if(!is.null(colnames(x@X))) colnames(x@X)[indices] else paste0("X", indices))
    } else if(type == "Y") {
        return(if(!is.null(colnames(x@Y))) colnames(x@Y)[indices] else paste0("Y", indices))
    }
}

#' Get sparsity information from sparse models
#' 
#' @param x A sparse model object
#' @param ... Additional arguments
#' @return List containing sparsity statistics
#' @export
sparsity_info <- function(x, ...) {
    UseMethod("sparsity_info")
}

setMethod("$", "SparseO2pls", function(x, name) {
    if(name %in% slotNames(x)) {
        return(slot(x, name))
    } else if(name %in% names(x@results)) {
        return(x@results[[name]])
    } else if(name %in% names(x@sparsity)) {
        return(x@sparsity[[name]])
    } else {
        stop("Unknown slot/element: ", name)
    }
})

setMethod("predict", "SparseO2pls", function(object, newdata, ...) {
    predict.SparseO2pls(object, newdata, ...)
})

#' @export
sparsity_info.sparse_o2pls <- function(x, ...) {
    return(list(
        selected_X = length(x$sparsity$selected_X),
        total_X = ncol(x$X),
        sparsity_X = x$sparsity$sparsity_X,
        percent_kept_X = round((1 - x$sparsity$sparsity_X) * 100, 1),
        selected_Y = length(x$sparsity$selected_Y),
        total_Y = ncol(x$Y),
        sparsity_Y = x$sparsity$sparsity_Y,
        percent_kept_Y = round((1 - x$sparsity$sparsity_Y) * 100, 1),
        keepX = x$sparsity$keepX,
        keepY = x$sparsity$keepY
    ))
}

#' @export
sparsity_info.sparse_plsda <- function(x, ...) {
    total_vars <- if(!is.null(x$total_variables)) x$total_variables else nrow(x$loadings)
    selected_vars <- length(x$selected_vars)
    sparsity <- 1 - selected_vars / total_vars
    
    return(list(
        selected = selected_vars,
        total = total_vars,
        sparsity = sparsity,
        percent_kept = round((1 - sparsity) * 100, 1),
        keepX = x$keepX
    ))
}

#' @export
sparsity_info.SparseO2pls <- function(x, ...) {
    return(list(
        selected_X = length(x@sparsity$selected_X),
        total_X = ncol(x@X),
        sparsity_X = x@sparsity$sparsity_X,
        percent_kept_X = round((1 - x@sparsity$sparsity_X) * 100, 1),
        selected_Y = length(x@sparsity$selected_Y),
        total_Y = ncol(x@Y),
        sparsity_Y = x@sparsity$sparsity_Y,
        percent_kept_Y = round((1 - x@sparsity$sparsity_Y) * 100, 1),
        keepX = x@sparsity$keepX,
        keepY = x@sparsity$keepY
    ))
}