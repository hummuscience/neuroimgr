##' Load data from suite2p
##'
##' Takes the output of a full suite2p run and turns it into a \linkS4class{CalciumExperiment} object
##' @title Import suite2p outputs
##' @param s2p.path The path to the folder containing suite2p outputs
##' @param time A time vector giving the timepoint at which every frame was taken. If left blank, a vector will be calculated with the fs parameter from Suite2p and the number of frames.
##' @param exp.n Experiment name/number. This is used to give the cells/ROIs unique IDs and will be randomly generated if left blank.
##' @param correction.factor The corection factor for neuropili substraction, defaults to 0.7
##' @param ... 
##' @return 
##' @author Muad Abd El Hay
##' @importFrom reticulate import
##' @export
readSuite2p <- function(s2p.path, time, exp.n, correction.factor=0.7, ...) {
  
  if(missing(exp.n)) {
  message("No experiment number/name given. Random name will be used.")
  exp.n <- randomAlphaNum()
  }

  np <- reticulate::import("numpy")
  os <- reticulate::import("os.path")
  
  message("Importing extracted fluorescence.")
  f <- np$load(
            os$expanduser(
                 paste(s2p.path,
                       "/plane0/F.npy",
                       sep = "")
               ),
            allow_pickle = TRUE)
  
  message("Importing neuropili traces.")
  fNeu <- np$load(
               os$expanduser(
                    paste(s2p.path,
                          "/plane0/Fneu.npy",
                          sep = "")
                  ),
               allow_pickle = TRUE) 

  isCell <- np$load(
                 os$expanduser(
                      paste(s2p.path,
                            "/plane0/iscell.npy",
                            sep = "")
                    ),
                 allow_pickle = TRUE)
  
  stats <- np$load(
                os$expanduser(
                     paste(s2p.path,
                           "/plane0/stat.npy",
                           sep = "")
                   ),
                allow_pickle = TRUE
              )

  message("Importing deconvoluted traces.")
  spks <- np$load(
               os$expanduser(
                    paste(s2p.path,
                          "/plane0/spks.npy",
                          sep = "")
                  ),
               allow_pickle = TRUE)

  ops <- np$load(
              os$expanduser(
                   paste(s2p.path,
                         "/plane0/ops.npy",
                         sep = "")
                 ),
              allow_pickle = TRUE)

  ops <- ops[[1]]
  
  message("Filtering cells.")
  fT <- t(f[isCell[,1]==1,])
  fNeuT <- t(fNeu[isCell[,1]==1,])
  spksT <- t(spks[isCell[,1]==1,])

  message(sprintf("Calculating corrected traces with factor:", correction.factor, sep=" "))
  fClean <- fT - correction.factor * fNeuT

  stats <- stats[which(isCell[,1] == 1)]

  message("Putting together cell parameters.")
  cellRadius <- purrr::map_dbl(stats, "radius")
  cellNpix <- purrr::map_dbl(stats, "npix")
  cellNpixNorm <- purrr::map_dbl(stats, "npix_norm")
  ## cellYpix <- purrr::map_dbl(stats, ~.x$ypix[1])
  ## cellXpix <- purrr::map_dbl(stats, ~.x$xpix[1])
  ## cellLam <- purrr::map_dbl(stats, ~.x$lam[1])
  cellCenterX <- purrr::map_dbl(stats, ~.x$med[[1]])
  cellCenterY <- purrr::map_dbl(stats, ~.x$med[[2]])
  cellCompact <- purrr::map_dbl(stats, "compact")
  cellFootprint <- purrr::map_dbl(stats, "footprint")
  cellAspectRatio <- purrr::map_dbl(stats, "aspect_ratio")
  cellSkew <- purrr::map_dbl(stats, "skew")
  cellStd <- purrr::map_dbl(stats, "std")
  
  phenoData <- DataFrame(radius = cellRadius,
                         npix = cellNpix,
                         npix_norm = cellNpixNorm,
                         ## ypix = cellYpix,
                         ## xpix = cellXpix,
                         ## lam = cellLam,
                         centerx <- cellCenterX,
                         centery <- cellCenterY,
                         compactness = cellCompact,
                         footprint = cellFootprint,
                         aspect_ratio = cellAspectRatio,
                         skewness = cellSkew,
                         std = cellStd)
  
  if (missing(time)) {
      message(sprintf("No time vector given. Generating time vector with framerate:",ops$fs,"Hz"))
      time <- 1:nrow(fT)/ops$fs
  } else {
    if (typeof(time) %in% c("integer","double")) {
      time = time
    } else {
      stop("Time vector supplied is neither double nor integer")  
    }     
  }

  ops <- c(exp.n, ops)

  message("Creating CalciumExperiment object.")
  ce <- CalciumExperiment(raw = fT,
                          colData = phenoData,
                          rowData = DataFrame(time),
                          metadata = ops)

  assay(ce, "neuropil") <- fNeuT
  assay(ce, "corrected") <- fClean
  assay(ce, "deconvoluted") <- spksT

  message("Generating unique cell names.")
  rownames(ce) <- paste("f", 1:nrow(ce), sep="")
  colnames(ce) <- paste(exp.n, 1:ncol(ce), sep="_")

  
  
  return(ce)
}
