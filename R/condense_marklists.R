#' Summarise cell lineage condition marker output in one csv file
#'
#' @param save_dir root directory used for running RNAseqEr. Default is getwd().
#' @param exclude_out_dir if there are folders in the out directory that should be
#' excluded, they can be specified here. An example is "data".
#' @param comp_mode Comparison mode. It cann be chosen between "overall" or
#' "clusterwise"
#' @param conditions Conditions that were tested such as Age, sex and tissue regiion
#' which may be profises as c("AgeGroup", "Sex", "Tissue)
#' @param cluster_label how is the cluster label supposed to be called in the output
#' file? The default is "cell_lineage".
#'
#' @return returns a list of marker genes where information of separate analyses are
#' collected in one data frame which is also saved to file
#' @export
#' @import dplyr
#'
#' @examples
#' Idents(cns) <- "AgeGroup"
#' Idents(oligodendrocytes) <- "AgeGroup"
#'
#' mark_cns <- FindAllMarkers(cns)
#' mark_ol <- FindAllMarkers(oligodendrocytes)
#'
#' save_dir_1 <- paste0(getwd(),
#'                             "/outs/all_celltypes/tables/condition_mark/cluster_id/overall")
#' save_dir_2 <- paste0(getwd(),
#'                             "/outs/oligodendrocytes/tables/condition_mark/cluster_id/overall")
#'
#' dir.create(save_dir_1, recursive = TRUE)
#' dir.create(save_dir_2, recursive = TRUE)
#'
#' write.csv(mark_cns, file.path(save_dir_1, "all_mark__AgeGroup_mark.csv")
#' write.csv(mark_ol, file.path(save_dir_2, "all_mark__AgeGroup_mark.csv")
#'
#' cond_marks <- condense_marklists(conditions = "AgeGroup")
#'
#'
condense_marklists <- function(save_dir = getwd(),
                               exclude_out_dir = "data",
                               comp_mode = "overall", #alternative is "clusterwise"
                               conditions,
                               cluster_label = "cell_lineage"){

  file_path = file.path(save_dir, "outs")

  dir_names_1 <- list.files(file_path)
  dir_name_exclude <- dir_names_1 %in% exclude_out_dir
  dir_names_1 <- dir_names_1[!dir_name_exclude]

  dir_names_2 <- list.files(file.path(file_path, dir_names_1, "tables",
                                      "condition_mark"))

  levels_dir_2 <- levels(as.factor(dir_names_2))

  for(l in 1:length(levels_dir_2)){
    fp_1 <-file.path(file_path, dir_names_1, "tables",
                     "condition_mark", levels_dir_2[l], comp_mode)
    bool <- dir.exists(fp_1)
    file_path_keep <- fp_1[bool]
    if(l == 1){
      file_path_seqEr <- file_path_keep
    }else{
      file_path_seqEr <- append(file_path_seqEr, file_path_keep)
    }
  }


  file_names <- list.files(file_path_seqEr)

  if(comp_mode ==  "overall"){
    fil_for_all <- grep("all", file_names)
    file_names <- file_names[fil_for_all]
  }

  for(l in 1:length(file_names)){
    fnames_1 <-file.path(file_path_seqEr, file_names[l])
    bool <- file.exists(fnames_1)
    file_keep <- fnames_1[bool]
    if(l == 1){
      file_path_use <- file_keep
    }else{
      file_path_use <- append(file_path_use, file_keep)

    }
  }


  csv_list <- list()


    for(j in 1: length(file_path_use)){
      if(file.exists(file_path_use[j]) == TRUE){
        curr_csv <- read.csv(file_path_use[j])

        curr_csv_name <- file_path_use[j]
        red_name <- strsplit(curr_csv_name, "/")
        red_name_2 <- red_name[[1]][length(red_name[[1]])]
        red_name_3 <- strsplit(red_name_2, ".csv")

        if(comp_mode ==  "clusterwise"){
        cond <- strsplit(red_name_3[[1]][1], "_")[[1]][1]
        cell_clu <- sub(paste0(cond, "_"), "", red_name_3)
        }

        if(comp_mode ==  "overall"){
          cond_split <- strsplit(red_name_3[[1]][1], "_")
          cond <- cond_split[[1]][length(cond_split[[1]])]

          red_na_pa_1 <- strsplit(curr_csv_name, "/outs/")[[1]][2]
          cell_clu <- strsplit(red_na_pa_1, "/tables/")[[1]][1]
        }

        curr_csv[["condition"]] <- cond
        curr_csv[["cluster_id"]] <- cell_clu

        csv_list[[j]] <- curr_csv

      }
    }

    df_to_save <- do.call("rbind", csv_list)

  save_dir_cre <- paste0(save_dir, "/outs/supplementary_tables/")
  if(dir.exists(save_dir_cre) == FALSE){
    dir.create(save_dir_cre, recursive = TRUE)
  }
  if(nrow(df_to_save) > 0){
    write.csv(df_to_save, file.path(save_dir_cre, paste0(comp_mode,
                                                    "_condition_markers.csv")))

  }

  return(df_to_save)
}
