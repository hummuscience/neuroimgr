##' Takes as input a CalciumExperiment object and a stimulus table that must contain the columns "stimulus", "begin", and "end" denoting the name/number of the stimulus, the beginning of the stimulus and the end, respectively. The stimulus table could either denote the times of stimuli using the units sipplied in the "time" field (check rowData of the CalciumExperiment object) or contain the frame numbers (row numbers) of the stimulus start and end. The default assumes a table contining the times.
##' @title Split CalciumExperiment object into multiple assays by stimuli
##' @param cexp CalciumExperiment to split.
##' @param stim.table The stimulus table with columns "stimulus", "begin", "end" where every row is one stimulus to split the dataset by.
##' @param stim.table.type The type of stimulus table that is supplied. 
##' @param buffer Whether to add a "buffer" before and after the stimulus. Useful for plotting.
##' @param buffer.size The relative size of the buffer. Defaults to half of the stimulus length.
##' @param ... 
##' @importFrom MultiAssayExperiment MultiAssayExperiment
##' @export
##' @return A MultiAssayExperiment object with each stimulus as independent CalciumExperiment object.
##' @author Muad Abd El Hay
splitByStimulus <- function(cexp, stim.table, stim.table.type = "time", buffer = FALSE, buffer.size = 0.5, ...){

  if (stim.table.type == "time") {
    
    frame.stim.table <- stim.table
    
    for (i in 1:nrow(frame.stim.table)){

      frame.stim.table[i, c("begin")] <- which.min(abs(rowData(cexp)$time - stim.table[i, c("begin")]))
      frame.stim.table[i, c("end")] <- which.min(abs(rowData(cexp)$time - stim.table[i, c("end")]))

    }
  } else if (stim.table.type == "frame") {
       
       frame.stim.table <- stim.table
       
  } else {

    stop("Invlid stim.table.type, must be either time or frame")

  }


  if (isTRUE(buffer)) {

    widths <- frame.stim.table$end - frame.stim.table$begin
    buffer.widths <- buffer.size * widths
    frame.stim.table$begin <- frame.stim.table$begin - buffer.widths
    frame.stim.table$end <- frame.stim.table$end + buffer.widths

    frame.stim.table$begin[which(frame.stim.table$begin < 0)] <- 0
    frame.stim.table$end[which(frame.stim.table$end > nrow(cexp))] <- nrow(cexp)

  }
  
  stimulusList <- vector(mode = "list", length = nrow(stim.table))
  names(stimulusList) <- frame.stim.table$stimulus
  
  for (i in 1:nrow(frame.stim.table)) {
    
    stimulusList[[i]] <- cexp[frame.stim.table[i,"begin"]:frame.stim.table[i,"end"],]

  }
  
  mcexp <- MultiAssayExperiment(experiments=stimulusList, 
                                colData = colData(stimulusList[[1]]))

  return(mcexp)
}

##' Function to estimate spikes and binarize a CalciumExperiment object.
##' @title Estimate spikes and binarize data using the L0 method by Jewell et al. 2019
##' @param cexp CalciumExperiment object to binarize.
##' @param slot The assay slot to use for binarization. Defaults to "normalized".
##' @param ... Additional parameters that can be passed to the spike estimation function.
##' @return
##' @export
##' @author Muad Abd El Hay
binarize <- function(cexp, slot = "normalized", ...){

  if (slot %in% names(assays(cexp))) {
    
    normalized.data <- assay(cexp, slot)

  } else {
    
    stop("no assay with given name found, please normalize the data first with the function normalize or assign another slot to perform binarization on.")
    
  }
  
  estimatedSpikes <- apply(normalized.data, 2, extractSpikesFromLZeroFit)

  colnames(estimatedSpikes) <- colnames(cexp)
  rownames(estimatedSpikes) <- rownames(cexp)
  
  assay(cexp, "l0spikes") <- estimatedSpikes

  return(cexp)
  
}

##' Utility function to extract the estimated spikes based on the FastLZeroSpikeInference package by Jewell et al. 2019.
##' @title Fit L0 spike inference model and extract spikes according to Jewell et al. 2019
##' @param x A vector containing the normalized fluorescence trace.
##' @param gam numeric. between 0 and 1. Parameter to tune spike inference.
##' @param lambda numeric. between 0 and 1. Parameter to tune spike inference.
##' @param constraint logical. Whether to run the fit with constrain or not.
##' @importFrom FastLZeroSpikeInference estimate_spikes
##' @return 
##' @author Muad Abd El Hay
extractSpikesFromLZeroFit <- function(x, gam = 0.8, lambda = 0.1, constraint = FALSE) {

  fit <- estimate_spikes(x, gam = gam, lambda = lambda, constraint = constraint)
  spikes <- rep(0, times = length(x))
  spikes[fit$spikes] <- 1

  return(spikes)
  
}

##' Function for fast center and scaling of matrix columns.
##' @title Fast center and/or scale columns using the matrixStats functions.
##' @param x 
##' @param center 
##' @param scale 
##' @param add_attr 
##' @param rows 
##' @param cols 
##' @return 
##' @references https://hopstat.wordpress.com/2016/02/23/a-faster-scale-function/
colScale = function(x,
    center = TRUE,
    scale = TRUE,
    add_attr = TRUE,
    rows = NULL,
    cols = NULL) {
 
    if (!is.null(rows) && !is.null(cols)) {
        x <- x[rows, cols, drop = FALSE]
    } else if (!is.null(rows)) {
        x <- x[rows, , drop = FALSE]
    } else if (!is.null(cols)) {
        x <- x[, cols, drop = FALSE]
    }
 
  ################
  # Get the column means
  ################
    cm = matrixStats::colMeans2(x, na.rm = TRUE)
  ################
  # Get the column sd
  ################
    if (scale) {
        csd = matrixStats::colSds(x, center = cm)
    } else {
        # just divide by 1 if not
        csd = rep(1, length = length(cm))
    }
    if (!center) {
        # just subtract 0
        cm = rep(0, length = length(cm))
    }
    x = t( (t(x) - cm) / csd )
    if (add_attr) {
        if (center) {
            attr(x, "scaled:center") <- cm
        }
        if (scale) {
            attr(x, "scaled:scale") <- csd
        }
    }
    return(x)
}

##' Function to perform exponential smoothing.
##' @title Exponential smoothing function.
##' @param x 
##' @param alpha 
##' @return 
##' @author Stephen Eglen
##' @references https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exposmoo.htm
expSmooth <- function(x, alpha=0.3) {
  y = x
  n = length(x)
  for (i in 2:n) {
    y[i] = (alpha*x[i]) + ((1-alpha)*y[i-1])
  }
  y
}

randomAlphaNum <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

##' Function to estimate baselines using the baseline package. Default method is "irls". For other options refer to the baseline function in the baseline package.
##' @title Function to estimate baselines of traces.
##' @param x 
##' @param method 
##' @param ... 
##' @return 
##' @author Muad Abd El Hay
##' @importFrom baseline baseline
estimateBaselines <- function(x, method = "irls", ...) {
  res <- baseline(t(x), method, ...)
  return(t(res@baseline))
}
