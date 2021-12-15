##' Plots staggered line plots of each ROI. Defaults to "raw" values without any additional modifications.
##'
##' Can be set to order the traces or even cluster them before plotting. Traces are scaled by default.
##' @title Function to plot CalciumExperiment objects.
#' @param x 
#' @param slot chracter. The assays() slot to use for plotting.
#' @param order logical. Whether to order the cells/ROIs by calculating the sums of the response.
#' @param cluster logical. Whether to clusted the cells/ROIs by correlation. 
##' @return A ggplot2 object with the traces for each cell.
##' @author Muad Abd El Hay
##' @export
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSums2
#' @import ggplot2
plotTraces <-  function(x, slot="raw", order=FALSE,cluster=FALSE, scale=TRUE) {
  if (slot == "raw"){
    chosen_assay <- raw(x)
  }
  else {
    chosen_assay <- SummarizedExperiment::assays(x)[[slot]]
  }
  
  if (isTRUE(scale)) {
    plotting_data <- colScale(chosen_assay)
  } else {
    plotting_data <- chosen_assay
  }

  if (isTRUE(cluster)) {
    cols.cor <- stats::cor(plotting_data, use = "pairwise.complete.obs", method = "pearson")
    cols.clust <- stats::hclust(stats::as.dist(1 - cols.cor))
    xorder <- cols.clust$order
    plotting_data <- plotting_data[,xorder]
  }

  if (isTRUE(order)) {
    xorder <- order(matrixStats::colSums2(plotting_data))
    plotting_data <- plotting_data[,xorder]
  }
  
  colnames(plotting_data) <- seq_len(ncol(plotting_data))
  
  color.factor <- ceiling(ncol(plotting_data)/9)

  plotting_data %>%
    tibble::as_tibble() %>%    
    dplyr::mutate(frame = seq_len(nrow(.))) %>%
    tidyr::gather(value = value, key = key, -frame) %>%
    dplyr::mutate(numkey = as.numeric(key)) %>% 
    ggplot(aes(x = frame, y = value + numkey, group = key, color = key)) +
    geom_line() +
    scale_colour_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),times=color.factor))+
    theme_void() +
    theme(legend.position = "none")
  
}

##' Plots staggered line plots of each ROI. Defaults to "raw" values without any additional modifications.
##'
##' Can be set to order the traces or even cluster them before plotting. Traces are scaled by default.
##' @title Function to plot CalciumExperiment objects.
#' @param x 
#' @param slot chracter. The assays() slot to use for plotting.
#' @param order logical. Whether to order the cells/ROIs by calculating the sums of the response.
#' @param cluster logical. Whether to clusted the cells/ROIs by correlation. 
##' @return A ggplot2 object with the traces for each cell.
##' @author Muad Abd El Hay
##' @export
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSums2
#' @import ggplot2
plotTraces2 <-  function(x, slot="raw", order=FALSE,cluster=FALSE, scale=TRUE) {
  if (slot == "raw"){
    chosen_assay <- raw(x)
  }
  else {
    chosen_assay <- SummarizedExperiment::assays(x)[[slot]]
  }
  
  if (isTRUE(scale)) {
    plotting_data <- colScale(chosen_assay)
  } else {
    plotting_data <- chosen_assay
  }

  if (isTRUE(cluster)) {
    cols.cor <- stats::cor(plotting_data, use = "pairwise.complete.obs", method = "pearson")
    cols.clust <- stats::hclust(stats::as.dist(1 - cols.cor))
    xorder <- cols.clust$order
    plotting_data <- plotting_data[,xorder]
  }

  if (isTRUE(order)) {
    xorder <- order(matrixStats::colSums2(plotting_data))
    plotting_data <- plotting_data[,xorder]
  }
  
  colnames(plotting_data) <- seq_len(ncol(plotting_data))
  
  color.factor <- ceiling(ncol(plotting_data)/9)

  plotting_data %>%
    tibble::as_tibble() %>%    
    dplyr::mutate(frame = seq_len(nrow(.))) %>%
    tidyr::gather(value = value, key = key, -frame) %>%
    dplyr::mutate(numkey = as.numeric(key)) %>% 
    ggplot(aes(x = frame, y = value + numkey, group = key)) +
    geom_path(aes(color = value)) +
    scale_colour_distiller(palette = "Spectral", direction = -1)+
    theme_void() +
    theme(legend.position = "none")
  
}

##' Plots staggered barcode plots of each ROI. Defaults to "l0spikes" values without any additional modifications.
##'
##' Can be set to order the traces or even cluster them before plotting. Traces are scaled by default.
##' @title Function to plot spikes from CalciumExperiment objects.
#' @param x 
#' @param slot chracter. The assays() slot to use for plotting.
#' @param order logical. Whether to order the cells/ROIs by calculating the sums of the response.
#' @param cluster logical. Whether to clusted the cells/ROIs by correlation. 
##' @return A ggplot2 object with the predicted pikes for each cell.
##' @author Muad Abd El Hay
##' @export
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSums2
#' @import ggplot2
plotSpikes <-  function(x, slot="l0spikes", order=FALSE,cluster=FALSE, scale=TRUE) {


  chosen_assay <- SummarizedExperiment::assays(x)[[slot]]
 
  
  if (isTRUE(scale)) {
    plotting_data <- colScale(chosen_assay)
  } else {
    plotting_data <- chosen_assay
  }

  if (isTRUE(cluster)) {
    cols.cor <- stats::cor(plotting_data, use = "pairwise.complete.obs", method = "pearson")
    cols.clust <- stats::hclust(stats::as.dist(1 - cols.cor))
    xorder <- cols.clust$order
    plotting_data <- plotting_data[,xorder]
  }

  if (isTRUE(order)) {
    xorder <- order(matrixStats::colSums2(plotting_data))
    plotting_data <- plotting_data[,xorder]
  }
  
  colnames(plotting_data) <- seq_len(ncol(plotting_data))
  
  color.factor <- ceiling(ncol(plotting_data)/9)

  plotting_data %>%
    tibble::as_tibble() %>%    
    dplyr::mutate(frame = seq_len(nrow(.))) %>%
    tidyr::gather(value = value, key = key, -frame) %>%
    dplyr::mutate(numkey = as.numeric(key)) %>% 
    dplyr::filter(value > 0) %>%
    ggplot(aes(x = frame, y = numkey, group = key, color = key)) +
    geom_point(shape = 108) +
    scale_colour_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),times=color.factor))+
    theme_void() +
    theme(legend.position = "none")
  
}

##' Plots staggered barcode plots of each ROI. Defaults to "l0spikes" values without any additional modifications.
##'
##' Can be set to order the traces or even cluster them before plotting. Traces are scaled by default.
##' @title Function to plot spikes from CalciumExperiment objects.
#' @param x 
#' @param slot chracter. The assays() slot to use for plotting.
#' @param order logical. Whether to order the cells/ROIs by calculating the sums of the response.
#' @param cluster logical. Whether to clusted the cells/ROIs by correlation. 
##' @return A ggplot2 object with the predicted pikes for each cell.
##' @author Muad Abd El Hay
##' @export
#' @importFrom magrittr %>%
#' @importFrom stats cor
#' @importFrom stats hclust
#' @importFrom stats as.dist
#' @importFrom rlang .data
#' @importFrom RColorBrewer brewer.pal
#' @importFrom matrixStats colSums2
#' @import ggplot2
plotSpikes2 <-  function(x, slot="l0spikes", order=FALSE,cluster=FALSE, scale=TRUE) {


  chosen_assay <- SummarizedExperiment::assays(x)[[slot]]
 
  
  if (isTRUE(scale)) {
    plotting_data <- colScale(chosen_assay)
  } else {
    plotting_data <- chosen_assay
  }

  if (isTRUE(cluster)) {
    cols.cor <- stats::cor(plotting_data, use = "pairwise.complete.obs", method = "pearson")
    cols.clust <- stats::hclust(stats::as.dist(1 - cols.cor))
    xorder <- cols.clust$order
    plotting_data <- plotting_data[,xorder]
  }

  if (isTRUE(order)) {
    xorder <- order(matrixStats::colSums2(plotting_data))
    plotting_data <- plotting_data[,xorder]
  }
  
  colnames(plotting_data) <- seq_len(ncol(plotting_data))
  
  color.factor <- ceiling(ncol(plotting_data)/9)

  plotting_data %>% 
    tibble::as_tibble() %>%    
    dplyr::mutate(frame = seq_len(nrow(.))) %>%
    tidyr::gather(value = value, key = key, -frame) %>%
    dplyr::mutate(numkey = as.numeric(key)) %>% 
    dplyr::filter(value > 0) %>%
    ggplot(aes(x = frame, y = numkey, group = key, color = key)) +
    geom_point(size = 0.1) +
    scale_colour_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),times=color.factor))+
    theme_void() +
    theme(legend.position = "none")
  
}
