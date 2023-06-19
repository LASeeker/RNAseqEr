#' Title
#' @description function performs DGE between all possible cluster combinations
#' at one or more given resolutions and saves a raw and filtered output
#' file for each comparison.
#' @param seur_obj quality controlled and clustered Seurat object
#' @param int_cols vector of column names that store cluster information
#' and will be used for all pairwise comparisons for clusters. This can be run
#' for one or more columns
#' @param only_pos TRUE/FALSE whether only positive log fold changes should be considered,
#' the default is TRUE.
#' @param min_pct minimal proportion of cells/nuclei that express a gene for it to
#' ne considered. Default is set higher here than in Seurat at 0.25.
#' @param logfc_threshold minimal log fold threshold between populations of
#' interest required for a gene to be considered for further testing.
#' @param fil_pct_1 minimum proportion of cells in the group of interest expressing
#' the cluster for filtered output
#' @param fil_pct_2 minimum proportion of cells in comparator group expressing
#' the cluster for filtered output
#' @param save_dir path to where output filed should be saved. default is working
#' directory
#' @param test_use statistical test to use for differential gene expression.
#' Default is "MAST" here but all other tests implemented in Seurat can be used
#' @param assay_use Default is RNA. In case integrated objects are used.
#'
#' @return saves lists of differentially expressed genes of pairwise compared
#' clusters. Save an unfiltered and a filtered list for each comparison
#' @import Seurat
#' @export
#'
#' @examples
#' library(Seurat)
#' pbmc_small_test <- seurat_proc(pbmc_small, tsne = FALSE)
#' pairwise_markers <- pairwise_dge(pbmc_small_test,
#'   int_cols = "RNA_snn_res.0.8"
#' )
#'
pairwise_dge <- function(seur_obj,
                         int_cols,
                         only_pos = TRUE,
                         min_pct = 0.25,
                         logfc_threshold = 0.25,
                         fil_pct_1 = 0.25,
                         fil_pct_2 = 0.1,
                         save_dir = getwd(),
                         test_use = "MAST",
                         assay_use = "RNA") {
  for (k in 1:length(int_cols)) {
    Idents(seur_obj) <- int_cols[k]
    # create all pairwise combinations at the picked resolution
    clust_id_list <- combn(levels(as.factor(seur_obj@meta.data[[int_cols]])), 2)
    for (i in 1:ncol(clust_id_list)) {
      clust_mark <- FindMarkers(seur_obj,
        ident.1 = clust_id_list[, i][[1]],
        ident.2 = clust_id_list[, i][[2]],
        min.pct = min_pct,
        test.use = test_use,
        assay = assay_use
      )
      clust_mark$cluster <- clust_id_list[, i][[1]]
      clust_mark$comp_to_clust <- clust_id_list[, i][[2]]
      write.csv(
        clust_mark,
        paste(save_dir,
          "/",
          int_cols[k],
          "_",
          clust_id_list[, i][[1]],
          "_",
          clust_id_list[, i][[2]],
          ".csv",
          sep = ""
        )
      )
    }
  }
}
