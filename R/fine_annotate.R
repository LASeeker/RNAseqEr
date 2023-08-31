#' Transfer fine clustering information from daughter datasets to mother dataset
#'
#' @param seur_obj Seurat object that shall be annotated
#' @param file_dir file directory where to find the processed daughter datasets.
#' Those have to include a column n the metadata that includes clustering information,
#' has the same name across all daughter datasets and is specified in the argument
#' cluster_col (see below).
#' @param cluster_col Metadata column name that is the same across all daughter
#' datasets and marks the column that includes cluster information that should
#' be added to the larger mother dataset (= seur_obj)
#' @param file_pattern This is the fie ending of the daughter datasets. The default
#' is ".RDS" but ".rds" will work too.
#'
#' @return returns a seurat object with added clustering information based on
#' subsets of the dataset.
#' @export
#' @import dplyr Seurat
#'
#' @examples
#' new_dir <- paste0(getwd(), "/example_data")
#' dir.create(new_dir)
#' saveRDS(oligodendrocytes, paste0(new_dir, "oligodendrocytes.RDS")
#' seur <- fine_annotate(cns,
#'                       file_dir = new_dir)
#'
fine_annotate <- function(seur_obj,
                          file_dir = getwd(),
                          cluster_col = "cluster_id",
                          file_pattern = ".RDS"
                          ){
  file_names <- list.files(file_dir, pattern = file_pattern)
  for(i in 1: length(file_names)){
    curr_file <- readRDS(paste(file_dir, file_names[i], sep = "/"))
    if(cluster_col %in% names(curr_file@meta.data)){
        new_df <- data.frame(Barcode = rownames(curr_file@meta.data),
                         RNAseqEr_cluster_id = curr_file@meta.data[[cluster_col]]
        )


        # if RNAseqER was used add potentially other interesting information
        RNAseqEr_cols <- grep("RNAseqEr", colnames(curr_file@meta.data))
        col_names <- colnames(curr_file@meta.data)[RNAseqEr_cols]
        if("RNAseqEr_annotation" %in% col_names){
          col_names<- col_names[!grepl(paste0("RNAseqEr_annotation", collapse = "|"), col_names)]
        }

        if(length(col_names) > 0){
          for(j in 1: length(col_names)){
            new_df[[col_names[j]]] <- curr_file@meta.data[[col_names[j]]]
          }

        }

        if(exists("save_df") == FALSE){
          save_df <- new_df
        }else{
          save_df <- rbind(save_df, new_df)
        }

    }
  }


  save_df_sort <- save_df[match(rownames(seur_obj@meta.data), save_df$Barcode),]

  for(y in 2: ncol(save_df_sort)){
    seur_obj@meta.data[[names(save_df_sort)[y]]] <- save_df_sort[,y]
  }

  DimPlot(seur_obj, group.by = "RNAseqEr_cluster_id", cols = colour_palette())
  DimPlot(seur_obj, group.by = "RNAseqEr_cluster_id", cols = colour_palette(),
          label = TRUE) + NoLegend()

  return(seur_obj)

}
