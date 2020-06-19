##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title 
#' @param x 
#' @param slot chracter. The assays() slot to use for plotting.
#' @param order logical. Whether to order the cells/ROIs by calculating the sums of the response.
#' @param cluster logical. Whether to clusted the cells/ROIs by correlation. 
##' @return 
##' @author Muad Abd El Hay
##' @export
#' @importFrom magrittr %>%
#' @import ggplot2
plotTraces <-  function(x, slot="raw", order=FALSE,cluster=FALSE) {
  if (slot == "raw"){
    chosen_assay <- raw(x)
  }
  else {
    chosen_assay <- assays(x)[[slot]]
  }
  
  scaled_data <- colScale(chosen_assay)
  
  if (isTRUE(cluster)) {
    cols.cor <- cor(scaled_data, use = "pairwise.complete.obs", method = "pearson")
    cols.clust <- hclust(as.dist(1 - cols.cor))
    xorder <- cols.clust$order
    scaled_data <- scaled_data[,xorder]
  }

  if (isTRUE(order)) {
    xorder <- order(matrixStats::colSums2(scaled_data))
    scaled_data <- scaled_data[,xorder]
  }
  
  colnames(scaled_data) <- 1:ncol(scaled_data)
  
  scaled_data %>%
    tibble::as_tibble() %>%    
    dplyr::mutate(frame = 1:nrow(.)) %>%
    tidyr::gather(value = value, key = key, -frame) %>%
    dplyr::mutate(numkey = as.numeric(key)) %>% 
    ggplot(aes(x = frame, y = value + numkey, group = key, color = key)) +
    geom_line() +
    theme_void() +
    theme(legend.position = "none")
  
}
