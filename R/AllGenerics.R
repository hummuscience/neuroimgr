##' Accessor function for raw fluorescence values per ROI
##' @title raw Accessor
#' @param x 
#' @param ... 
##' @return 
##' @author Muad Abd El Hay
#' @export
setGeneric("raw", function(x, ...) standardGeneric("raw"))


##' Accessor function for raw fluorescence values per ROI
##' @title raw Accessor
#' @param x 
#' @param withDimnames 
##' @return 
##' @author Muad Abd El Hay
#' @export
#' @importFrom SummarizedExperiment assay
setMethod("raw", "CalciumExperiment", function(x, withDimnames=TRUE) {
    assay(x, withDimnames=withDimnames)
})

#' @title Function to normalize traces.
#' @export
#' @importFrom BiocGenerics normalize
##' @param object CalciumExperiment object. 
##' @param method character. The method to use for baseline correction. Defaults to "estimateBaseline",
##' @param slot character. Which assay slot to use for baseline correction.
##' @param window numeric. The width of the baseline window to use for dF/F0 calculation. A window of width 500 will take the first 500 frames as baseline.
##' @param ... 
setMethod("normalize", "CalciumExperiment", function(object, method = "estimateBaseline", slot, window, ...) {
  
  non.normalized.data <- assay(object,slot)
  if (missing(slot)) {
    stop("Chose which assay to normalize by setting the slot argument.")
  }

  if (method == "estimateBaseline") {
    baseline <- estimateBaselines(non.normalized.data)
    dff <- (non.normalized.data - baseline) / baseline
  }

  if (method == "start") {
    if (missing(window)) {
      stop("No baseline starting window specified.")
    }
    baseline <- apply(non.normalized.data[1:window,],2, mean)
    df <- sweep(non.normalized.data, 2, baseline, "-")
    dff <- sweep(df, 2, baseline, "/")
  }

  assay(object, "normalized") <- dff

  return(object)

}
)
