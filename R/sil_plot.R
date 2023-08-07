#' Creating silhouette plots for different clustering resolutions
#'
#' @param seur_obj Quality controlles, if necessary batch corrected and integrated
#' Seurat Object that has been clustered at different resolutions
#' @param reduction dimensional reduction that should be used. Default is "PCA"
#' @param col_pattern A specific pattern of column names that contain clustering
#' information. Default is "RNA_snn_res"
#' @param plot_cols specify colours for plots. Default is colour_palette()
#' @param clust_lab TRUE/FALSE if should clisters be labeled
#' @param save_to_file TRUE/FALSE whether silhouette plots should be saved to file
#' in specified save_dir or only plotted in environment
#' @param save_dir specify where output should be saved. Defaut is working directory
#' @param dir_lab label used for which data is analysed. Default is "all_celltypes"
#' @param width width of output plot, default is 7
#' @param height height of output plot, default is 5
#'
#' @return creates visual output of silhouettes to be inspected for cluster
#' stability which are saved to file in save_dir and returns a dataframe with
#' all silhouette data of different resolutions combined
#' @import Seurat SingleCellExperiment gtools cluster
#' @export
#'
#' @examples
#' library(Seurat)
#' library(SingleCellExperiment)
#' library(gtools)
#' library(cluster)
#' cns <- seurat_proc(cns,tsne = FALSE)
#' av_sil_df <- sil_plot(cns)
sil_plot <- function(seur_obj,
                     reduction = "PCA",
                     col_pattern = "RNA_snn_res",
                     plot_cols = colour_palette(),
                     clust_lab = TRUE,
                     save_to_file = FALSE,
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
  distance <- dist(reducedDim(sce_obj, reduction))
  if(save_to_file == TRUE){
    sil_dir <- paste0(save_dir, "/outs/", dir_lab, "/plots/")
  }


  for(i in 1: length(names_col)){
    clust <- met_dat[[names_col[i]]]

    if(length(levels(clust)) > 1){
            clust_int <- as.integer(paste0(clust))
            sil <- silhouette(clust_int, dist = distance)
            sil_cols <- plot_cols[ifelse(sil[,3] > 0, sil[,1]+1, sil[,2]+1)]
            sil_cols <- sil_cols[order(-sil[,1], sil[,3])]

            if(save_to_file == TRUE){
            pdf(paste0(sil_dir,
                       names_col[i], "_sil.pdf"), width=width, height=height)
            plot(sil, border = NA)
            plot(sil,
                 main = paste("clusters"),
                 border = sil_cols,
                 col = sil_cols,
                 do.col.sort = FALSE)
            dev.off()
            }

            plot(sil,
                 main = paste("clusters"),
                 border = sil_cols,
                 col = sil_cols,
                 do.col.sort = FALSE)

          if(exists("av_sil_df") == FALSE){
            av_sil_df <- data.frame(res = names_col[i],
                                    av_sil_w = summary(sil)$avg.width)
          }else{
            append_df <- data.frame(res = names_col[i],
                                    av_sil_w = summary(sil)$avg.width)
            av_sil_df <- rbind(av_sil_df, append_df)
            }

            if(length(levels(clust)) == 1 | length(levels(clust)) == 0){
                print("More than one cluster required. Increase cluster resolution")
            }
    }
  }
  return(av_sil_df)
}

