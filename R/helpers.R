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
    
    for (i in seq_len(nrow(frame.stim.table))){

      frame.stim.table[i, c("begin")] <- which.min(abs(SummarizedExperiment::rowData(cexp)$time - stim.table[i, c("begin")]))
      frame.stim.table[i, c("end")] <- which.min(abs(SummarizedExperiment::rowData(cexp)$time - stim.table[i, c("end")]))

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
  
  for (i in seq_len(nrow(frame.stim.table))) {
    
    stimulusList[[i]] <- cexp[seq.int(from = frame.stim.table[i, "begin"], to = frame.stim.table[i,"end"]) ,]

  }
  
  mcexp <- MultiAssayExperiment(experiments=stimulusList, 
                                colData = SummarizedExperiment::colData(stimulusList[[1]]))

  return(mcexp)
}

##' Utility function to extract the estimated spikes based on the FastLZeroSpikeInference package by Jewell et al. 2019.
##' @title Fit L0 spike inference model and extract spikes according to Jewell et al. 2019
##' @param x A vector containing the normalized fluorescence trace.
##' @param gam numeric. between 0 and 1. Parameter to tune spike inference.
##' @param lambda numeric. between 0 and 1. Parameter to tune spike inference.
##' @param constraint logical. Whether to run the fit with constrain or not.
##' @importFrom FastLZeroSpikeInference estimate_spikes
##' @return A vector containig the estimated spikes as 1s and 0s.
##' @author Muad Abd El Hay
extractSpikesFromLZeroFit <- function(x, gam = 0.8, lambda = 0.1, constraint = FALSE) {

  fit <- FastLZeroSpikeInference::estimate_spikes(x, gam = gam, lambda = lambda, constraint = constraint)
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
##' @return A matrix with scaled columns.
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

######
#' Calculate the Area Under Curve of y~x
#'
#'@param y Your y values (measures ?)
#'@param x Your x values (time ?)
#'@param start : The first x value 
#'@param stop : The last x value
#'@param na.stop : returns NA if one value is NA
#'@param ex.na.stop : returns NA if the first or the last value is NA
## getAUC <- function(y, x, start=S4Vectors::head(x,1), stop=S4Vectors::tail(x,1), na.stop=FALSE, ex.na.stop=TRUE){
##   if(all(is.na(y))) return(NA)
##   bounds = seq_along(which(x==start):which(x==stop))
##   x=x[bounds]
##   y=y[bounds]
##   r = which(is.na(y))
##   if(length(r)>0){
##     if(na.stop==TRUE) return(NA)
##     if(ex.na.stop==TRUE & (is.na(S4Vectors::first(y)) | is.na(last(y)))) return(NA)
##     if(is.na(last(y))) warning("Last value is NA, so this AUC is bad and you should feel bad", call. = FALSE) 
##     if(is.na(S4Vectors::first(y))) warning("First value is NA, so this AUC is bad and you should feel bad", call. = FALSE) 
##     x = x[-r]
##     y = y[-r]
##   }
##   sum(diff(x[order(x)])*zoo::rollmean(y[order(x)],2))
## }

##' Function to perform exponential smoothing.
##' @title Exponential smoothing function.
##' @param x 
##' @param alpha 
##' @return 
##' @author Stephen Eglen
##' @references https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/exposmoo.htm
## expSmooth <- function(x, alpha=0.3) {
##   y = x
##   n = length(x)
##   for (i in seq_along(2:n)) {
##     y[i] = (alpha*x[i]) + ((1-alpha)*y[i-1])
##   }
##   y
## }

##' Helper function to calculate a firing rate estimate using a Bayesian Adaptive Kernel Smoother as suggested by Ahmadi et al. 2018
##' https://doi.org/10.1371/journal.pone.0206794
##' @title Firing Rate Estimation using Bayesian Adaptive Kernel Smoother
##' @param spikes A binarized vector (0 and 1) of spikes
##' @param t A vector of times in seconds, same length as spikes
##' @param a shape parmeter according to Ahmadi et al., defaults to 4
##' @param b scale parameter according to Ahmadi et al., defaults to 0.975
##' @return A vector of firing rate density estimation
##' @author Muad Abd El Hay
baks <- function(spikes, t, a = 4, b = 0.975) {

  if( sum(spikes) == 0) {

    firing_rate = rep(0, length(t))

  } else {

    spike_times <- t[spikes == 1]
    
    n <- length(spike_times)


    
    b <- b**n

    sumnum <- 0

    sumdenom <-  0

    for (i in seq_len(n)) {

      numerator = (((t - spike_times[i])**2)/2 + 1/b)**(-a)
      denominator = (((t - spike_times[i])**2)/2 +1/b)**(-a - 0.5)
      sumnum = sumnum + numerator
      sumdenom = sumdenom + denominator
      h = (gamma(a)/gamma(a + 0.5)) * (sumnum/sumdenom)

      firing_rate = rep(0, length(t))
    }

    for (j in seq_len(n)) {
      k = (1 / (sqrt(2 * pi) * h)) * exp(-((t - spike_times[j])**2) / (2*h**2))
      firing_rate = firing_rate + k                           
    }
  }
  return(firing_rate)
  
}

##' Helper function to create an assay name when none is supplied.
##' @title Create a random alpha-numeric string of length n.
##' @param n Length of the alphanumeric string.
##' @return A string with random alpha-numerics of length n.
##' @author Muad Abd El Hay
randomAlphaNum <- function(n = 1) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

##' Function to estimate baselines using the baseline package. Default method is "irls". For other options refer to the baseline function in the baseline package.
##' @title Function to estimate baselines of traces.
##' @param x 
##' @param method 
##' @param ... 
##' @return A baseline-corrected trace.
##' @author Muad Abd El Hay
##' @importFrom baseline baseline
estimateBaselines <- function(x, method = "modpolyfit", ...) {
  res <- baseline::baseline(t(x), method, ...)
  return(t(res@baseline))
}
