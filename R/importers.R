##' Load data from suite2p
##'
##' Takes the output of a full suite2p run and turns it into a \linkS4class{CalciumExperiment} object
##' @title Import suite2p outputs
##' @param s2p.path The path to the folder containing suite2p outputs
##' @param time A time vector giving the timepoint at which every frame was taken. If left blank, a vector will be calculated with the fs parameter from Suite2p and the number of frames.
##' @param exp.n Experiment name/number. This is used to give the cells/ROIs unique IDs and will be randomly generated if left blank.
##' @param correction.factor The corection factor for neuropili substraction, defaults to 0.7
##' @param ... 
##' @return A CalciumExperiment object with the data from a Suite2p output folder.
##' @author Muad Abd El Hay
##' @importFrom reticulate import
##' @export
readSuite2p <- function(s2p.path, time, exp.n, correction.factor=0.7, filter = TRUE, ...) {
  
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
  
  if (isTRUE(filter)) {

  message("Filtering cells.")
  fT <- t(f[isCell[,1]==1,])
  fNeuT <- t(fNeu[isCell[,1]==1,])
  spksT <- t(spks[isCell[,1]==1,])

    } else {
  message("Taking all ROIs.")    
  fT <- t(f)
  fNeuT <- t(fNeu)
  spksT <- t(spks)    

}

  
    
  if (missing(time)) {
    message(sprintf("No time vector given. Generating time vector with framerate:",ops$fs,"Hz"))
    time <- seq_len(nrow(fT))/ops$fs
  } else {
    if (typeof(time) %in% c("integer","double")) {
      time = time
    } else {
      stop("Time vector supplied is neither double nor integer")  
    }     
  }

  if(!(nrow(fT) == length(time))) {

    stop("Time vector length is not qual to number of frames.")
    
  }
  

  message(sprintf("Calculating corrected traces with factor:", correction.factor, sep=" "))
  fClean <- fT - correction.factor * fNeuT

  if (isTRUE(filter)) {

  stats <- stats[which(isCell[,1] == 1)]
}
  phenoData <-  S4Vectors::DataFrame(cellid = paste(exp.n, seq_len(ncol(fClean)), sep="_"))

  message("Putting together cell parameters.")
  if("radius" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, radius = purrr::map_dbl(stats, "radius"))}
  if("npix" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, npix = purrr::map_dbl(stats, "npix"))}
  if("npix_norm" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, npix_norm = purrr::map_dbl(stats, "npix_norm"))}
  if("med" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, centerx = purrr::map_dbl(stats, ~.x$med[[1]]))
    phenoData <- cbind(phenoData, centery = purrr::map_dbl(stats, ~.x$med[[2]]))}
  if("compact" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, compactness = purrr::map_dbl(stats, "compact"))}
  if("footprint" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, footprint = purrr::map_dbl(stats, "footprint"))}
  if("aspect_ratio" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, aspect_ratio = purrr::map_dbl(stats, "aspect_ratio"))}
  if("skew" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, skewness = purrr::map_dbl(stats, "skew"))}
  if("std" %in% names(stats[[1]])) {
    phenoData <- cbind(phenoData, std = purrr::map_dbl(stats, "std"))}
  
  xpix <-  purrr::map(stats, "xpix")
  ypix <-  purrr::map(stats, "ypix")
  xypix <- mapply(cbind, xpix, ypix, SIMPLIFY=F)
  xypix <- lapply(xypix, as.data.frame)

  xypix <- Map('+', xypix, 1)
  
  tidymasks <- dplyr::bind_rows(setNames(xypix, seq_along(xypix)), .id = "id")
  
  masks <- matrix(0, nrow = ops$Lx, ncol = ops$Ly)
  
  for (i in seq_len(length(xypix))) {

    cell.xy.pix <- xypix[[i]]
    
    masks[as.matrix(cell.xy.pix)] <- i
  }

  ops$exp.n <- exp.n
  ops$tidymasks <- tidymasks
  ops$masks <- masks
  

  message("Creating CalciumExperiment object.")
  ce <- CalciumExperiment(raw = fT,
                          colData = phenoData,
                          rowData = S4Vectors::DataFrame(time),
                          metadata = ops)

  SummarizedExperiment::assay(ce, "neuropil") <- fNeuT
  SummarizedExperiment::assay(ce, "corrected") <- fClean
  SummarizedExperiment::assay(ce, "deconvoluted") <- spksT

  message("Generating unique cell names.")
  rownames(ce) <- paste("f", seq_len(nrow(ce)), sep="")
  colnames(ce) <- paste(exp.n, seq_len(ncol(ce)), sep="_")

  
  
  return(ce)
}
