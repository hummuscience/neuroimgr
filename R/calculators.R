##' Function that takes a MultiAssayExperiment object containings CalciumExperiment objects and calculates parameters for each assay (stimulus).
##' @title Calculate stimulus parameters from a split CalciumExperiment (i.e. MultiAssayExperiment) object.
##' @param mcexp MultiAssayExperiment object containing multiple CalciumExperiment objects (one per stimulus) to calculate from.
##' @param ... 
##' @return A MultiAssayExperiment object with calculated base parameters for each CalciumExperiment object.
##' @export
##' @author Muad Abd El Hay
calculateParameters <- function(mcexp, ...) {

  for (i in seq_len(length(mcexp))) {

    mcexp[[i]] <- getBaseParameters(mcexp[[i]])
  }

return(mcexp)

}

##' Calculates base parameters of a CalciumExperiment object such as mean, max, min responses as well as number of predicted spikes and location of first spike. 
##' @title Calculate base parameters. 
##' @param cexp CalciumExperiment object to calculate base parameters from. Should be normalized and binarized before running.
##' @param slot 
##' @param binslot 
##' @param ... 
##' @return A CalciumExperiment object with calculated base parameters added to the colData DataFrame
##' @author Muad Abd El Hay
getBaseParameters <-  function(cexp, slot = "normalized", binslot = "l0spikes", ...) {

  if(!(is(cexp, "CalciumExperiment"))) {

    stop("cexp must be of class CalciumExperiment")

  }
  
  if (!("normalized" %in% names(SummarizedExperiment::assays(cexp)))) {
    warning("Chosen assay for parameter calculation is not normalized, are you sure you want to continue?")
  }

  if (!(binslot %in% names(SummarizedExperiment::assays(cexp)))) {
    stop("Chosen binarized assay for parameter calculation is not present, forgot to binarize?")
  }

  normalized.traces <- SummarizedExperiment::assay(cexp, slot)
  binarized.traces <- SummarizedExperiment::assay(cexp, binslot)
  
  stim_mean <- apply(normalized.traces, 2, mean)
  stim_max <- apply(normalized.traces, 2, max)
  stim_min <- apply(normalized.traces, 2, min)
  spikenum <- apply(binarized.traces,2, sum)
  first_spike <- rep(NA, times = ncol(cexp))

  for (i in seq_len(ncol(cexp))) {

    first_spike[i] <- Position( function(x) x > 0, binarized.traces[,i])
  }
  
  base.parameters <- S4Vectors::DataFrame(stim_mean,
                               stim_max,
                               stim_min,
                               spikenum,
                               first_spike)
  
  SummarizedExperiment::colData(cexp) <- cbind(SummarizedExperiment::colData(cexp), base.parameters)
  
  return(cexp)

}

##' Calculates stimulus-dependent parameters of a CalciumExperiment object such as stimulus threshold, correlation, and mutual information.
##' @title Calculate base parameters. 
##' @param cexp CalciumExperiment object that has base parameters calculated already. Should be normalized and binarized before running.
##' @param slot The assay slot representing the normalized calcium traces
##' @param binslot The assay slot with the binarized traces
##' @param stimulus.trace The rowData column containing the stimulus to calculate for
##' @param ... 
##' @return A CalciumExperiment object with calculated base parameters added to the colData DataFrame
##' @author Muad Abd El Hay
##' @export
getStimulusParameters <-  function(cexp, slot = "normalized", binslot = "l0spikes", stimulus.trace, ...) {

  if(!(is(cexp, "CalciumExperiment"))) {

    stop("cexp must be of class CalciumExperiment")

  }
  
  if (!("normalized" %in% names(SummarizedExperiment::assays(cexp)))) {
    warning("Chosen assay for parameter calculation is not normalized, are you sure you want to continue?")
  }

  if (!(binslot %in% names(SummarizedExperiment::assays(cexp)))) {
    stop("Chosen binarized assay for parameter calculation is not present, forgot to binarize?")
  }

  normalized.traces <- SummarizedExperiment::assay(cexp, slot)
  binarized.traces <- SummarizedExperiment::assay(cexp, binslot)
  stimulus <-  SummarizedExperiment::rowData(cexp)[,stimulus.trace]

  stimulus.threshold <- rep(NA, times = ncol(cexp))

  for (i in seq_len(ncol(cexp))) {
    if (is.na(SummarizedExperiment::colData(cexp)$first_spike[i])) next
    stimulus.threshold[i] <- stimulus[SummarizedExperiment::colData(cexp)$first_spike[i]]
  }
  
  
  stimulus.parameters <- S4Vectors::DataFrame(stimulus.threshold)
  
  SummarizedExperiment::colData(cexp) <- cbind(SummarizedExperiment::colData(cexp), stimulus.parameters)
  
  return(cexp)

}

##' Function to estimate spikes and binarize a CalciumExperiment object.
##' @title Estimate spikes and binarize data using the L0 method by Jewell et al. 2019
##' @param cexp CalciumExperiment object to binarize.
##' @param slot The assay slot to use for binarization. Defaults to "normalized".
##' @param ... Additional parameters that can be passed to the spike estimation function.
##' @return A CalciumExperiment object with an additional assay containig the estimated spikes as 0 and 1 binary.
##' @export
##' @author Muad Abd El Hay
binarize <- function(cexp, slot = "normalized", ...){

  if (slot %in% names(SummarizedExperiment::assays(cexp))) {
    
    normalized.data <- SummarizedExperiment::assay(cexp, slot)

  } else {
    
    stop("no assay with given name found, please normalize the data first with the function normalize or assign another slot to perform binarization on.")
    
  }
  
  estimatedSpikes <- apply(normalized.data, 2, extractSpikesFromLZeroFit)

  colnames(estimatedSpikes) <- colnames(cexp)
  rownames(estimatedSpikes) <- rownames(cexp)
  
  SummarizedExperiment::assay(cexp, "l0spikes") <- estimatedSpikes

  return(cexp)
  
}

##' A wrapper function for firing rate estimation of binarized CalciumExperiment objects.
##' @title Firing rate estimation from binarized data
##' @param cexp CalciumExperiment object that underwent binarization
##' @param slot Assay slot containing binarized data
##' @param ... Variables passed to the bayesian adaptive kernel smoother (see \link{baks} for more information).
##' @return A CalciumExperiment object with firing rate estimation
##' @author Muad Abd El Hay
##' @export
estimateFiringRate <- function(cexp, slot = "l0spikes", ...){

  if (slot %in% names(SummarizedExperiment::assays(cexp))) {
    
    binarized.data <- SummarizedExperiment::assay(cexp, slot)

  } else {
    
    stop("no assay with given name found, please binarize the data first with the function binarize or assign another slot to perform firing rate estimation on.")
    
  }
  
  estimatedFiringRate <- apply(binarized.data, 2, baks, t = rowData(cexp)$time)

  colnames(estimatedFiringRate) <- colnames(cexp)
  rownames(estimatedFiringRate) <- rownames(cexp)
  
  SummarizedExperiment::assay(cexp, "frs") <- estimatedFiringRate

  return(cexp)
  
}
