#' Generate heatmaps
#' @description
#' Caculate average gene expression of a provided list of genes across different
#' clusters and plot in heatmaps that are saved to file.
#'
#' @param seur_obj Quality controlled and if necessary batch corrected and
#' integrated Seurat object with several clustering resolutions
#' @param use_resol TRUE/FALSE whether column labels should be used that indicate
#' standard names for different clustering resolutions opposed to user indicated
#' cluster labels that for example specify a cell type. Default is TRUE.
#' @param col_names if use_resol = FALSE a string of at least one column names
#' has to be provided.
#' @param col_pattern Which column pattern to use to find clustering resolutions.
#' Default is"RNA_snn_res"
#' @param int_genes list of interesting genes that should be tested.
#' @param gene_count if list of interesting genes is long it will be shortened
#' to gene_count. Default is 250.
#' @param save_dir where output should be saved, default is working directory
#' @param label TRUE/FALSE wheter cluster labels should be displayed on plot.
#' Default is TRUE
#' @param plot_cols colours for plot, default is colour_palette from the RNAseqEr
#' package
#' @param draw_lines TRUE/FALSE whether white lines should be drawn between columns
#' corresponding to each cluster. Default is FALSE
#' @param width width of output plot. Default is 5.
#' @param height heigth of output plot. Default is 15.
#' @param size size of labels above bar. Default is 11.
#' @param diff_threshold resolution threshold for determining if clusters are
#' different enough. The average expression is considered.
#'
#' @return saves plots to file
#' @export
#'
#' @examples
#' library(Seurat)
#' seur <- seurat_proc(pbmc_small,tsne = FALSE)
#' int_genes <- head(VariableFeatures(object = seur), 250)
#' heatmap_seqEr(seur, int_genes = int_genes)
heatmap_seqEr <- function(seur_obj,
                          use_resol = TRUE,
                          col_names,
                          col_pattern = "RNA_snn_res",
                          int_genes,
                          gene_count = 200,
                          save_dir = getwd(),
                          label = TRUE,
                          plot_cols = colour_palette(),
                          draw_lines = FALSE,
                          width=3000,
                          height = 10000,
                          x_size = 11,
                          y_size = 15,
                          diff_threshold = 10){
  if(use_resol == TRUE){
    res_col <- grep(pattern = col_pattern, names(seur_obj@meta.data))
    names_col <- names(seur_obj@meta.data)[res_col]
    # gtools function, sorts gene_names alphanumeric:
    names_col <- mixedsort(names_col)
  }else{
    names_col <- col_names
  }


  met_dat <- as.data.frame(seur_obj@meta.data)

  if(length(int_genes > gene_count)){
    int_genes_sh <- paste(head(int_genes, gene_count))
  }else{
    int_genes_sh <- int_genes
  }


  for(i in 1: length(names_col)){
    if(length(levels(as.factor(met_dat[[names_col[i]]])))>1){


    cluster_averages <- AverageExpression(seur_obj,
                                          group.by = names_col[i],
                                          return.seurat = TRUE)
    cluster_averages@meta.data$cluster <- levels(as.factor(met_dat[[names_col[i]]]))
    cluster_averages@meta.data$cluster <- factor(cluster_averages@meta.data$cluster,
                                                 levels = levels(as.factor(met_dat[[names_col[i]]])))


    hm_av <- DoHeatmap(object = cluster_averages,
                       features = int_genes_sh,
                       label = label,
                       group.by = "cluster",
                       group.colors = plot_cols,
                       draw.lines = draw_lines,
                       size = x_size)

    hm_av <- hm_av + theme_minimal(y_size)

    png(paste0(save_dir,
               "/",
               names_col[i],
               "_average_marker_heatmap.png"),
        width=width, height=height, res = 300)

    print(hm_av)

    dev.off()

    print(hm_av)


    # compare custers in object
    cluster_averages_df <- as.data.frame(AverageExpression(seur_obj,
                                                           group.by = names_col[i],
                                                           return.seurat = FALSE))
    colnames(cluster_averages_df) <- levels(as.factor(met_dat[[names_col[i]]]))

    clust_id_mtx <- combn(levels(as.factor(met_dat[[names_col[i]]])), 2)

    for(k in 1: ncol(clust_id_mtx)){
      ident_1 = clust_id_mtx[, k][[1]]
      ident_2 = clust_id_mtx[, k][[2]]

      max_diff <- abs(max(cluster_averages_df[ident_1] - cluster_averages_df[ident_2]))
      if(max_diff < diff_threshold){
        Print(paste0("When considering resolution/ clustering ",
                     names_col[i],
                     ", ",
                     clust_id_mtx[, k][[1]],
                     "and ",
                     ident_2 = clust_id_mtx[, k][[2]],
                     " should be merged as they appaear too similar."))
      }else{
        print(max_diff)
      }

    }


    }
  }
  print("Generated heatmap and saved them to file.")
}
