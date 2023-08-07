#' plot_list
#' @description Function that plots dimensionally reduced data at different
#' resolutions and saves plot as outputs
#' @param seur_obj Seurat object with dimensional reductions and several resolutions.
#' @param col_pattern name pattern used to identify columns that contain clustering
#' information at different resolutions. Default is "RNA_snn_res."
#' @param plot_cols color palette default is colour_pallette() from the RNAseqEr
#' package
#' @param clust_lab Whether DimPlot should be annotated with clusters labels
#' @param label_size size of cluster labels
#' @param save_dir file path to where plots should be saved, default is working
#' directory
#' @param dir_lab Label indicating if all celltypes are being processed (default =
#' "all_celltypes) or a cell lineage which is important for output folder structure.
#' @param width width of pdf file
#' @param height height of pdf file
#' @param use_reduction which dimensional reduction should be used for plotting
#' ("pca". "tsne", "umap)
#' @import Seurat gtools
#' @return function creates output folders where it saves DimPlots at different
#' resolutions as pdf
#' @export
#'
#' @examples
#' library(Seurat)
#' library(RColorBrewer)
#' library(ggsci)
#' cns <- seurat_proc(cns, tsne = FALSE)
#' plot_list(cns)
plot_list <- function(seur_obj,
                      col_pattern = "RNA_snn_res.",
                      plot_cols = colour_palette(),
                      clust_lab = TRUE,
                      label_size = 8,
                      save_dir = getwd(),
                      dir_lab = "all_celltypes",
                      width = 7,
                      height = 5,
                      use_reduction = "umap") {
  extr_res_col <- grep(pattern = col_pattern, names(seur_obj@meta.data))
  res_names <- names(seur_obj@meta.data[extr_res_col])
  # gtools function, sorts gene_names alphanumeric:
  res_names <- gtools::mixedsort(res_names)
  save_dir_path <- (paste0(
      save_dir,
      "/outs/",
      dir_lab,
      "/plots/resolution_plots/"))
  dir.create(save_dir_path, recursive = TRUE)

  plot_l <- list()
  for (i in 1:length(res_names)) {
    pdf(paste0(save_dir_path,
      res_names[i], "_umap.pdf"
    ), width = width, height = height)
    dim_plot <- DimPlot(seur_obj,
      reduction = use_reduction,
      cols = plot_cols,
      group.by = res_names[i],
      label = clust_lab,
      label.size = label_size
    ) + NoLegend()
    print(dim_plot)
    dev.off()
    print(dim_plot)
  }
}
