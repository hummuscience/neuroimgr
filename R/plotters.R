##' Plots staggered line plots of each ROI. Defaults to "raw" values without any additional modifications.
##'
##' Can be set to order the traces or even cluster them before plotting. Traces are scaled by default.
##' @title Function to plot CalciumExperiment objects.
#' @param x 
#' @param slot chracter. The assays() slot to use for plotting.
#' @param order logical. Whether to order the cells/ROIs by calculating the sums of the response.
#' @param cluster logical. Whether to clusted the cells/ROIs by correlation. 
##' @return 
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
    chosen_assay <- assays(x)[[slot]]
  }
  
  if (isTRUE(scale)) {
    plotting_data <- colScale(chosen_assay)
  } else {
    plotting_data <- chosen_assay
  }

  if (isTRUE(cluster)) {
    cols.cor <- cor(plotting_data, use = "pairwise.complete.obs", method = "pearson")
    cols.clust <- hclust(as.dist(1 - cols.cor))
    xorder <- cols.clust$order
    plotting_data <- plotting_data[,xorder]
  }

  if (isTRUE(order)) {
    xorder <- order(matrixStats::colSums2(plotting_data))
    plotting_data <- plotting_data[,xorder]
  }
  
  colnames(plotting_data) <- 1:ncol(plotting_data)
  
  color.factor <- ceiling(ncol(plotting_data)/9)

  plotting_data %>%
    tibble::as_tibble() %>%    
    dplyr::mutate(frame = 1:nrow(.)) %>%
    tidyr::gather(value = value, key = key, -frame) %>%
    dplyr::mutate(numkey = as.numeric(key)) %>% 
    ggplot(aes(x = frame, y = value + numkey, group = key, color = key)) +
    geom_line() +
    scale_colour_manual(values=rep(RColorBrewer::brewer.pal(9,"Set1"),times=color.factor))+
    theme_void() +
    theme(legend.position = "none")
  
}
