#' int_res_all_mark
#' @description
#' This function performs a overall differential gene expression analysis using
#' Seurat's FindAllMarkers() at a chosen number of resoltions and saves both
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
#' @param assay_use assay to use for differential gene expression analysis.
#' Default is "RNA". For sketched data, use "sketch".
#' @param dir_lab specifies which cell lineage is investigated which is important for
#' output folder structure. The default is "all_celltypes".
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
                             test_use = "MAST",
                             assay_use = "RNA",
                             dir_lab = "all_celltypes",
                             min_cells_per_cluster = 3,
                             min_cells_for_comparison = 5,
                             verbose = TRUE) {
  for (i in 1:length(int_cols)) {
    curr_res <- int_cols[i]
    
    # Validate cluster sizes before analysis
    validation <- validate_cluster_sizes(seur_obj, curr_res, 
                                        min_cells_per_cluster, min_cells_for_comparison)
    
    if(verbose) {
      cat("=== Overall DGE Analysis for", curr_res, "===\n")
      cat("Total cells:", validation$total_cells, "\n")
      cat("Number of clusters:", validation$n_clusters, "\n")
      cat("Valid clusters:", validation$n_valid_clusters, "\n")
    }
    
    # Check if analysis is feasible
    if(!validation$feasible) {
      if(verbose) {
        cat("❌ Skipping overall DGE for", curr_res, ": insufficient cells/clusters\n")
      }
      next
    }
    
    # Filter small clusters if needed
    if(validation$recommendations$action == "filter_clusters") {
      if(verbose) {
        cat("⚠️  Filtering small clusters for", curr_res, "\n")
      }
      filtered_obj <- filter_small_clusters(seur_obj, curr_res, min_cells_per_cluster, verbose)
      if(is.null(filtered_obj)) {
        if(verbose) {
          cat("❌ Skipping overall DGE for", curr_res, ": no valid clusters after filtering\n")
        }
        next
      }
      analysis_obj <- filtered_obj
    } else {
      analysis_obj <- seur_obj
    }
    
    Idents(analysis_obj) <- curr_res
    
    if(length(levels(Idents(analysis_obj))) > 1){
      tryCatch({
        all_mark <- FindAllMarkers(analysis_obj,
          only.pos = only_pos,
          min.pct = min_pct,
          logfc.threshold = logfc_threshold,
          test.use = test_use,
          assay = assay_use,
          verbose = FALSE
        )
        fil_mark <- subset(
          all_mark,
          all_mark$pct.1 > fil_pct_1 &
            all_mark$pct.2 < fil_pct_2
        )

        save_path <- paste0(save_dir,
                            "/outs/",
                            dir_lab,
                            "/tables/cluster_marker/overall"
                            )

        if(dir.exists(save_path) == FALSE){
          dir.create(save_path, recursive = TRUE)
        }

        write.csv(all_mark, paste(save_path, "/all_mark", curr_res, ".csv", sep = ""))
        write.csv(fil_mark, paste(save_path, "/fil_mark", curr_res, ".csv", sep = ""))
        
        if(verbose) {
          cat("✅ Overall DGE completed for", curr_res, "\n")
          cat("   Markers found:", nrow(all_mark), "\n")
          cat("   Filtered markers:", nrow(fil_mark), "\n")
        }
        
      }, error = function(e) {
        if(verbose) {
          cat("❌ Error in overall DGE for", curr_res, ":", e$message, "\n")
        }
      })
    } else {
      if(verbose) {
        cat("⚠️  Skipping overall DGE for", curr_res, ": only one cluster\n")
      }
    }
  }
}
