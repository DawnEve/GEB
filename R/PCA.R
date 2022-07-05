

#' Plot cell PC score as color on UMAP plot
#'
#' @param sce Seurat object
#' @param nPCs an integer, total number of PC to plot;
#'
#' @return a ggplot2 merge obj, by cowplot::plot_grid
#' @export
#' @importFrom purrr map
#' @importFrom Seurat FetchData
#' @import ggplot2 magrittr cowplot
#'
#' @examples
#' plot_PC_on_UMAP(pbmc, ncol=4)
#' plot_PC_on_UMAP(pbmc, 9)
plot_PC_on_UMAP=function(sce, nPCs=16, ncol=3){
  # Defining the information in the seurat object of interest
  columns <- c( paste0("PC_", 1:nPCs),  "seurat_clusters", "UMAP_1", "UMAP_2")

  # Extracting this data from the seurat object
  pc_data <- FetchData(sce, vars = columns)

  # Center of cluster on UMAP
  umap_label_pos=FetchData(sce, vars = c("seurat_clusters", "UMAP_1", "UMAP_2")) %>%
    group_by(seurat_clusters) %>%
    summarise(x=mean(UMAP_1), y=mean(UMAP_2))

  # Plotting a UMAP plot for each of the PCs
  purrr::map(paste0("PC_", 1:nPCs), function(pc){
    ggplot(pc_data, aes(UMAP_1, UMAP_2)) +
      geom_point(aes_string(color=pc), alpha = 0.7) +
      #scale_color_gradient(guide = FALSE, low = "grey90", high = "red")  +
      scale_color_gradient2(low = "blue", mid="white", high = "red")  +
      geom_text(data=umap_label_pos, aes(label=seurat_clusters, x, y)) +
      ggtitle(pc)+ theme_bw()
  }) %>%
    cowplot::plot_grid(plotlist = ., ncol=ncol)
}
