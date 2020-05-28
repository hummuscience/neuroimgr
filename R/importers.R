##' Load data from suite2p
##'
##' Takes the output of a full suite2p run and turns it into a \linkS4class{CalciumExperiment} object
##' @title Import suite2p outputs
##' @param s2p.path The path to the folder containing suite2p outputs
##' @param correction.factor The corection factor for neuropili substraction, defaults to 0.7
##' @return 
##' @author Muad Abd El Hay
##' @importFrom reticulate import
readSuite2p <- function(s2p.path, correction.factor=0.7) {

  np <- reticulate::import("numpy")
  os <- reticulate::import("os.path")
  
  f <- np$load(
            os$expanduser(
                 paste(s2p.path,
                       "/plane0/F.npy",
                       sep = "")
               ),
            allow_pickle = TRUE)
  
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

  
  spks <- np$load(
               os$expanduser(
                    paste(s2p.path,
                          "/plane0/spks.npy",
                          sep = "")
                  ),
               allow_pickle = TRUE)

  f <- f[isCell[,1]==1,]
  fNeu <- fNeu[isCell[,1]==1,]
  spks <- spks[isCell[,1]==1,]
  fClean <- f - correction.factor * fNeu

  stats <- stats[which(isCell[,1] == 1)]

  cellRadius <- purrr::map_dbl(stats, "radius")
  cellNpix <- purrr::map_dbl(stats, "npix")
  cellNpixNorm <- purrr::map_dbl(stats, "npix_norm")
  cellYpix <- purrr::map_dbl(stats, ~.x$ypix[1])
  cellXpix <- purrr::map_dbl(stats, ~.x$xpix[1])
  cellLam <- purrr::map_dbl(stats, ~.x$lam[1])
  cellCompact <- purrr::map_dbl(stats, "compact")
  cellFootprint <- purrr::map_dbl(stats, "footprint")
  cellAspectRatio <- purrr::map_dbl(stats, "aspect_ratio")
  cellSkew <- purrr::map_dbl(stats, "skew")
  cellStd <- purrr::map_dbl(stats, "std")
  
  phenoData <- DataFrame(radius = cellRadius,
                         npix = cellNpix,
                         npix_norm = cellNpixNorm,
                         ypix = cellYpix,
                         xpix = cellXpix,
                         lam = cellLam,
                         compactness = cellCompact,
                         footprint = cellFootprint,
                         aspect_ratio = cellAspectRatio,
                         skewness = cellSkew,
                         std = cellStd)
  
  ce <- CalciumExperiment(t(f),
                          colData = phenoData)

  assay(ce, "neuropil") <- t(fNeu)
  assay(ce, "corrected") <- t(fClean)
  assay(ce, "deconvoluted") <- t(spks)

  
  
  return(ce)
}
