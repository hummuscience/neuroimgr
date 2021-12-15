##' Function that defines the CalciumExperiment class
##'
##' 
##' @title CalciumExperiment class
#' @param raw 
#' @param ... 
##' @return 


##' @author Muad Abd El Hay
##'@export
#' @import methods
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
.CalciumExperiment <- setClass("CalciumExperiment", contains = "SummarizedExperiment")

##' Constructor for CalciumExperiment objects.
##'
##' 
##' @title CalciumExperiment constructor
#' @param raw 
#' @param ... 
##' @return A CalciumExperiment object.
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
