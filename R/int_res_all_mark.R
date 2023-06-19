#' int_res_all_mark
#' @description
#' This function performs a overall differential gene expression analysis using
#' Seurat's FindAllMarkers() at a chosen number of resultions and saves both
#' raw and filtered lists in a provided output directory.
#' @param seur_obj Seurat object that has been quality controlled and clustered.
#' @param int_cols list of colums of interest that should be compared with one another
#' using differential gene expression analyses.
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
#' @return saves lists of differentially expressed genes between clusters of
#' interest to output folders once filtered once unfiltered.
#' @import Seurat MAST
#' @export
#'
#' @examples
#' library(Seurat)
#' pbmc_small_test <- seurat_proc(pbmc_small, tsne = FALSE)
#' int_res_all_mark(pbmc_small_test, int_cols = c(
#'   "RNA_snn_res.0.8",
#'   "RNA_snn_res.1"
#' ))
int_res_all_mark <- function(seur_obj,
                             int_cols,
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
        test.use = test_use
      )
      fil_mark <- subset(
        all_mark,
        all_mark$pct.1 > fil_pct_1 &
          all_mark$pct.2 < fil_pct_2
      )


      write.csv(all_mark, paste(save_dir, "/all_mark", int_cols[i], ".csv", sep = ""))
      write.csv(fil_mark, paste(save_dir, "/fil_mark", int_cols[i], ".csv", sep = ""))
    }
  }
}
