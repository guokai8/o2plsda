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