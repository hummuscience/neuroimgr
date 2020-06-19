##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
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
