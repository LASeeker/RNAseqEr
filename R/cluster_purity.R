#' Calculating cluster purity measurements and plotting them to the environment
#' and file
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
#' @param dir_lab label used for which data is analysed. Default is "all_celltypes"
#' @param width width of output plot, default is 7
#' @param height height of output plot, default is 5
#'
#' @return saves plots to files.
#' @export
#' @import Seurat SingleCellExperiment gtools cluster bluster ggplot2 ggbeeswarm
#'
#' @examples
#' library(Seurat)
#' library(SingleCellExperiment)
#' library(gtools)
#' library(cluster)
#' library(bluster)
#' library(ggplot2)
#' library(ggbeeswarm)
#' cns <- seurat_proc(cns,tsne = FALSE)
#' clu_pure(cns)
clu_pure <- function(seur_obj,
                     reduction = "PCA",
                     col_pattern = "RNA_snn_res",
                     plot_cols = colour_palette(),
                     clust_lab = TRUE,
                     label_size = 8,
                     save_dir = getwd(),
                     dir_lab = "all_celltypes",
                     width=7,
                     height=5){
  sce_obj <- as.SingleCellExperiment(seur_obj)
  res_col <- grep(pattern = col_pattern, names(colData(sce_obj)))
  names_col <- names(colData(sce_obj))[res_col]
  # gtools function, sorts gene_names alphanumeric:
  names_col <- mixedsort(names_col)
  met_dat <- as.data.frame(colData(sce_obj))
  clu_pure_dir <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/tables/cluster_purity_data/")
  dir.create(clu_pure_dir, recursive = TRUE)

  clu_pure_plot_dir <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/plots/cluster_purity_plots/")
  dir.create(clu_pure_plot_dir, recursive = TRUE)

  for(i in 1: length(names_col)){
    clust <- met_dat[[names_col[i]]]
    if(length(levels(clust)) > 1){
    clust_int <- as.integer(paste0(clust))

    pure <- neighborPurity(reducedDim(sce_obj, reduction), clusters = clust_int)
    pure_data <- as.data.frame(pure)
    pure_data$maximum <- factor(pure_data$maximum)
    pure_data$cluster <- factor(clust_int)

    write.csv(pure_data, paste0(clu_pure_dir,
                                "/",
                                names_col[i],
                                "_cluster.purity.csv"))


    pure_plot <- ggplot(pure_data, aes(x=cluster, y=purity, colour=maximum)) +
      ggbeeswarm::geom_quasirandom(method="smiley") +
      theme_bw(20) +
      xlab(names_col[i])

    pdf(paste0(clu_pure_plot_dir,
               names_col[i],
               "_clu_purity.pdf"),
        width=width,
        height=height)

    print(pure_plot)

    dev.off()

    print(pure_plot)
    }
  }
  print("Done with calculating and plotting cluster purity measures at selected resolutions")
}
