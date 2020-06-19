##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
#' @param raw 
#' @param ... 
##' @return 
##' @author Muad Abd El Hay
##' #' @export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.CalciumExperiment <- setClass("CalciumExperiment", contains = "SummarizedExperiment")

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
#' @param raw 
#' @param ... 
##' @return 
##' @author Muad Abd El Hay
##' @export
#' @importFrom SummarizedExperiment SummarizedExperiment
CalciumExperiment <- function(raw, ...) {
  se <- SummarizedExperiment(list(raw=raw), 
                             ...)
  .CalciumExperiment(se)
}

## setValidity2("CalciumExperiment", function(object) {
##   msg <- NULL
##   if (is.null(msg)) {
##     TRUE
##   } else msg
## })
