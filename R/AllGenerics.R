##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
#' @param x 
#' @param ... 
##' @return 
##' @author Muad Abd El Hay
#' @export
setGeneric("raw", function(x, ...) standardGeneric("raw"))


##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
#' @param x 
#' @param withDimnames 
##' @return 
##' @author Muad Abd El Hay
#' @export
#' @importFrom SummarizedExperiment assay
setMethod("raw", "CalciumExperiment", function(x, withDimnames=TRUE) {
    assay(x, withDimnames=withDimnames)
})
