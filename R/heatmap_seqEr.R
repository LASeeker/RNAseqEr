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
#' @param x_size size of labels above bar. Default is 11.
#' @param y_size Legend size. Default is 15.
#' @param max_diff_threshold resolution threshold for determining if clusters are
#' different enough based on the maximum difference of average gene expression between
#' two clusters.Default is 10 (differences must be larger than threshold).
#' @param mean_diff_thres resolution threshold for determining if clusters are
#' different enough based on the mean difference of average gene expression between
#' two clusters.Default is 0.01 (differences must be larger than threshold).
#' @param eucl_dist_thres Threshold of Eucledian distance difference between
#' the average expression across all int_genes when comparing two clusters to determine
#' if they are different enough or too similar and should be merged. Default is
#' 30.
#' @param print_values TRUE/FALSE whether max gene expression difference, mean
#' gene expresiion difference and Eucleadian discance should be printed for the
#' comparison of each two cluster pair that passes the above set thresholds. Defautl
#' is FALSE. Setting it to TRUE may enable finding meaningful thresholds.
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
                          max_diff_threshold = 10,
                          mean_diff_thres = 0.1,
                          eucl_dist_thres = 30,
                          print_values = FALSE){
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
    y = 0
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


    # compare clusters in object
    cluster_averages_df <- as.data.frame(AverageExpression(seur_obj,
                                                           group.by = names_col[i],
                                                           return.seurat = FALSE))
    colnames(cluster_averages_df) <- levels(as.factor(met_dat[[names_col[i]]]))

    clust_id_mtx <- combn(levels(as.factor(met_dat[[names_col[i]]])), 2)


    for(k in 1: ncol(clust_id_mtx)){
      ident_1 = clust_id_mtx[, k][[1]]
      ident_2 = clust_id_mtx[, k][[2]]

      all_diff <- abs(cluster_averages_df[ident_1] - cluster_averages_df[ident_2])
      names(all_diff) <- "diff"

      max_diff <- abs(max(cluster_averages_df[ident_1] - cluster_averages_df[ident_2]))

      euclidean <- function(a, b) sqrt(sum((a - b)^2))
      eucl_dist <- euclidean(cluster_averages_df[ident_1], cluster_averages_df[ident_2])



      if(max_diff < max_diff_threshold &
         mean(all_diff$diff, na.rm = TRUE) < mean_diff_thres |
         eucl_dist < eucl_dist_thres){
        print(paste0("When considering resolution/ clustering ",
                     names_col[i],
                     ", ",
                     clust_id_mtx[, k][[1]],
                     "and ",
                     ident_2 = clust_id_mtx[, k][[2]],
                     " may need to be merged as they may be too similar based on set thresholds."))
        y = y + 1
      }else{
        if(print_values == TRUE){
          print(c("When considering resolution/ clustering",
                     names_col[i],
                     ", and the comparison of ",
                     ident_1,
                     "and",
                     ident_2,
                     "the max difference is ",
                     max_diff,
                     ", the mean difference is ",
                     mean(all_diff$diff, na.rm = TRUE),
                     "and the Eucledian distance is ",
                     eucl_dist,
                     "."
                    ))
        }


      }

    }


    }
    df <- data.frame(resolution = names_col[i],
                     cluster_sim_score = y)

    if(i == 1){
      keep_df <- df
    }else{
      keep_df <- rbind(keep_df, df)
    }
  }
  print("Generated heatmap and saved them to file.")
  return(keep_df)
}
