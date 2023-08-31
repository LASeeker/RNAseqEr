#' Subset Seurat object based on cluster QC data
#' @description
#' This function assumes that RNAseqEr's cluster_qc() was previously run on
#' the data. Alternatively, it works with other metadata columns that contain
#' quality control metrics, as long as the clusters to keep are indicated in the
#' "filter" argument. If RNAseqEr is used (RNASeqEr = TRUE), file paths will be
#' constructed based on previous steps using save_dir as root folder. If other
#' QC data is used, set RNASeqEr to FALSE
#' @param seur_obj Seurat object
#' @param cluster_col metadata column that contains clustering information. Default is
#' "RNAseqEr_annotation"
#' @param save_dir root directory used for analysis.
#' @param dir_lab label of cell lineage. Default is "all_celltypes"
#' @param filter Filter indicating which clusters should be retained when
#' considering the "conclusion" columns in the QC output. The default is "Pass"
#' which means that all other clusters are going to be removed.
#' @return Seurat object that has some clusters removed baesed on previous
#' cluster quality control.
#' @export
#' @examples
#' here()
#' seur_obj <- cluster_qc(cns,
#'                    cluster_col = "rough_annot",
#'                    RNAseqEr = TRUE,
#'                    save_dir = here(),
#'                    vars = c("caseNO", "process_number"))
#'
#' seur_obj <- remove_clu(seur_obj,
#'                        cluster_col = "rough_annot",
#'                        save_dir = here())
remove_clu <- function(seur_obj,
                       cluster_col = "RNAseqEr_annotation",
                       RNAseqEr = TRUE,
                       save_dir = getwd(),
                       dir_lab = "all_celltypes",
                       filter = "Pass"){
  if(RNAseqEr == TRUE){
  file_path <- paste0(save_dir, "/outs/", dir_lab, "/tables/cluster_qc")
  } else if(RNAseqEr == FALSE){
    file_path <- save_dir()
  }else{
    print("Set RNAseqEr to TRUE if cluster_qc() has been used for cluster_qc() and
          to false if alternative cluster QC approaches were used.")
  }
  files <- list.files(file_path)
  my_files <- lapply(paste(file_path, files, sep = "/"), read.csv)
  my_df <- as.data.frame(my_files)
  col_numb <- grep("conclusion", names(my_df))

  col_names <- names(my_df[col_numb])

  for(i in 1:length(col_names)){
    bool <- my_df[[col_names[i]]] != filter
    subs_df <- my_df[bool, ]
    cluster <- subs_df[cluster_col]

    if(i == 1){
      save_df <- cluster
    }else{
      save_df <- rbind(save_df, cluster)
    }
  }

  uniq_cluster <- unique(save_df)

  Idents(seur_obj) <- cluster_col

  seur_obj <-

  test <- subset(seur_obj, idents = uniq_cluster[[cluster_col]], invert = TRUE)

  DimPlot(seur_obj, group.by = cluster_col)

  print(paste0(uniq_cluster[[cluster_col]], " cluster was removed because it failed set QC thresholds."))

  return(seur_obj)




}
