#' Aproximate Silhouette
#'
#' @param seur_obj Quality controlled and if necessary batch corrected and
#' integrated Seurat object with several clustering resolutions
#' @param reduction which dimensional reduction to use. Default is "PCA"
#' @param col_pattern Which column pattern to use to find clustering resolutions.
#' Default is"RNA_snn_res"
#' @param plot_cols colours for plot, default is colour_palette from the RNAseqEr
#' package
#' @param clust_lab TRUE/FALSE whether cluster labels should be displayed
#' @param label_size size of cluster labels
#' @param save_dir where output should be saved, default is working directory
#' @param width width of output plot, default is 7
#' @param height height of output plot, default is 5
#'
#' @return calculates approximate silhouette and plots results to file
#' @import Seurat SingleCellExperiment gtools cluster bluster ggplot2 ggbeeswarm
#' @export
#'
#' @examples
#' library(Seurat)
#' library(SingleCellExperiment)
#' library(gtools)
#' library(cluster)
#' library(bluster)
#' library(ggplot2)
#' library(ggbeeswarm)
#' seur <- seurat_proc(pbmc_small,tsne = FALSE)
#' approx_sil(seur)
#'
approx_sil <- function(seur_obj,
                       reduction = "PCA",
                       col_pattern = "RNA_snn_res",
                       plot_cols = colour_palette(),
                       clust_lab = TRUE,
                       label_size = 8,
                       save_dir = getwd(),
                       width=7,
                       height=5){
  sce_obj <- as.SingleCellExperiment(seur_obj)
  res_col <- grep(pattern = col_pattern, names(colData(sce_obj)))
  names_col <- names(colData(sce_obj))[res_col]
  # gtools function, sorts gene_names alphanumeric:
  names_col <- mixedsort(names_col)
  met_dat <- as.data.frame(colData(sce_obj))

  for(i in 1: length(names_col)){
    clust <- met_dat[[names_col[i]]]
    if(length(levels(clust)) > 1){
      clust_int <- as.integer(paste0(clust))

      sil_approx <- approxSilhouette(reducedDim(sce_obj, reduction),
                                     clusters = clust_int)
      sil_data <- as.data.frame(sil_approx)
      sil_data$closest <- factor(ifelse(sil_data$width > 0, clust_int, sil_data$other))
      sil_data$cluster <- factor(clust_int)

      apr_sil_plot <-ggplot(sil_data, aes(x=cluster, y=width, colour=closest)) +
        ggbeeswarm::geom_quasirandom(method="smiley") + theme_bw(20) +
        xlab(names_col[i])


      pdf(paste0(save_dir, "/",
                 names_col[i], "_sil.pdf"), width=width, height=height)


      print(apr_sil_plot)

      dev.off()

      print(apr_sil_plot)
    }
  }
  print("Done with approximate silhouette calculation and plotting")
}
