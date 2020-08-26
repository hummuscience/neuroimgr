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

##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
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
