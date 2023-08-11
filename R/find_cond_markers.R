#' Find condition markers
#' @description
#' By using differential gene expression analysis over the entire dataset and
#' within clusters, markers for all conditions of interest (disease/control,
#' treatment/control, young/old, tissue regions,...) are identified and saved.
#' @param seur_obj Seurat object that has been quality controlled and clustered.
#' @param cluster_id metadata column that contains cluster information used.
#' @param int_cols list of metadata columns that contain all conditions of interest
#' (such as disease/control, treatment/control, young/old, tissue regions,...).
#' If condition of interest is a factor with more than 2 levels, pairwise
#' comparisons are implemented.
#' @param dir_lab label to organise output into cell linages.
#' Default is "all_celltypes".
#' @param only_pos whether only positive log-fold changes should be considered,
#' default is TRUE.
#' @param min_pct minimum proportion within cluster of interest expressing the gene.
#' @param logfc_threshold threshold for log2 fold changes
#' @param fil_pct_1 filtering threshold for the minimum proportion of cells within
#' the cluster of interest expressing the gene.
#' @param fil_pct_2 filtering threshold for the maximum proportion of cells outside
#' the cluster of interest expressing the gene.
#' @param save_dir path to folder where output files should be saved, default is
#' current working directory
#' @param test_use statistical approach being used for differential gene expression
#' analysis. Default is MAST.
#'
#' @return saves lists of differentially expressed genes between conditions to
#' folders
#' @import Seurat MAST
#' @export
#'
#' @examples
#' library(Seurat)
#' cns <- seurat_proc(cns, tsne = FALSE)
#' find_cond_markers(cns,
#'                   int_cols = c("AgeGroup","Tissue"),
#'                   cluster_id = "rough_annot")
find_cond_markers <- function(seur_obj,
                             cluster_id,
                             int_cols,
                             dir_lab = "all_celltypes",
                             only_pos = TRUE,
                             min_pct = 0.25,
                             logfc_threshold = 0.25,
                             fil_pct_1 = 0.25,
                             fil_pct_2 = 0.6,
                             save_dir = getwd(),
                             test_use = "MAST") {
  for (i in 1:length(int_cols)) {
    Idents(seur_obj) <- int_cols[i]
    if(length(levels(Idents(seur_obj))) > 1){
      all_mark <- FindAllMarkers(seur_obj,
                                 only.pos = only_pos,
                                 min.pct = min_pct,
                                 logfc.threshold = logfc_threshold,
                                 test.use = test_use,
                                 verbose = FALSE
      )
      fil_mark <- subset(
        all_mark,
        all_mark$pct.1 > fil_pct_1 &
          all_mark$pct.2 < fil_pct_2
      )

      dir_save <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/tables/condition_mark/",
                         cluster_id,
                         "/overall")

      if(dir.exists(dir_save) == FALSE){
        dir.create(dir_save, recursive = TRUE)
      }


      write.csv(all_mark, paste(dir_save, "/all_mark__", int_cols[i], ".csv", sep = ""))
      write.csv(fil_mark, paste(dir_save, "/fil_mark__", int_cols[i], ".csv", sep = ""))
    }
    for(k in 1: length(levels(as.factor(seur_obj@meta.data[[cluster_id]])))){
      curr_clu <- levels(as.factor(seur_obj@meta.data[[cluster_id]]))[k]
      Idents(seur_obj) <- cluster_id
      curr_clu_srt <- subset(seur_obj, ident = curr_clu)

      Idents(curr_clu_srt) <- int_cols[i]
      within_clu_mark <- FindAllMarkers(curr_clu_srt,
                                 only.pos = only_pos,
                                 min.pct = min_pct,
                                 logfc.threshold = logfc_threshold,
                                 test.use = test_use,
                                 verbose = FALSE)


      dir_save_clu <- paste0(save_dir,
                         "/outs/",
                         dir_lab,
                         "/tables/condition_mark/",
                         cluster_id,
                         "/clusterwise")

      if(dir.exists(dir_save_clu) == FALSE){
        dir.create(dir_save_clu, recursive = TRUE)
      }

      if(nrow(within_clu_mark > 0)){

          write.csv(within_clu_mark, paste0(dir_save_clu,
                                            "/",
                                            int_cols[i],
                                            "_",
                                            curr_clu,
                                            ".csv"))
        }

    }

  }
}
