##' Class "O2pls"
##' This class represents the Annotation information
##' @name O2pls-class
##' @aliases O2pls-class
##' @docType class
##' @slot X a Numeric matrix (input)
##' @slot Y a Numeric matrix (input)
##' @slot params paramaters ysed in o2pls analysis
##' @slot results list of o2pls results
##' @exportClass O2pls
##' @author Kai Guo
##' @keywords classes
setClass("O2pls",
         representation = representation(
             X = "matrix",
             Y = "matrix",
             params ="list",
             results ="list"
         ))


##' Class "SparseO2pls"
##' This class represents sparse O2PLS analysis results
##' @name SparseO2pls-class
##' @aliases SparseO2pls-class
##' @docType class
##' @slot X Numeric matrix (input X)
##' @slot Y Numeric matrix (input Y)
##' @slot params List of parameters including sparsity
##' @slot results List of sparse O2PLS results
##' @slot sparsity List of sparsity-related information
##' @slot preprocessing list of preprocessing information
##' @slot call The matched call used to generate the object
##' @exportClass SparseO2pls
##' @author Kai Guo
setClass("SparseO2pls",
         representation = representation(
             X = "matrix",
             Y = "matrix", 
             params = "list",
             results = "list",
             sparsity = "list",
             preprocessing = "list",
             call = "call"
         ))


##' Class "TuneResult"
##' @exportClass TuneResult
setClass("TuneResult",
         representation = representation(
             optimal = "list",
             all_results = "list", 
             method = "character",
             measure = "character",
             validation = "character"
         ))

setClass("sparse_plsda",
         representation = representation(
             call = "call",
             ncomp = "numeric",
             keepX = "numeric",
             scores = "matrix",
             loadings = "matrix",
             vip = "matrix",
             selected_vars = "numeric",
             classification = "list",
             preprocessing = "list",
             dist = "character",
             Y_levels = "character",
             X = "matrix",
             Y = "factor"
         ))

setClass("stability_selection",
         representation = representation(
             selection_probabilities = "array",
             stable_variables = "list",
             threshold = "numeric",
             n_bootstrap = "numeric",
             keepX_range = "numeric",
             method = "character",
             parameters = "list"
         ))
