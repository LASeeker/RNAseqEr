#' max_pure
#' @description
#' This function uses cluster purity measurements at different resoution and looks
#' for a resolution suitable for a rough celltype annotattion.
#'
#' @param save_path path to where cluster purity measures are saved as
#' .csv files. Default is "getwd()/outs/all_celltypes/tables/cluster_purity_data".
#' @param exponent this is used to weigh the cluster purity measure. A quotient of
#' cluster purity ^ exponent / cluster count is used to calculate an overall pure
#' measure.
#' @param pure_thres Threshold of cluster purity to be used. The function is looking
#' for the maximum number of clusters over this purity threshold.
#'
#' @return
#' A recommendation of a clustering resolution that maximizes cluster purity and
#' cluster number and saves summary table to new folder within read_dir.
#' @export
#'
#' @examples
#' library(here)
#' save_dir_cr <- here("outs", "example_out")
#' dir.create(save_dir_cr, recursive = TRUE)
#'
#' cns <- seurat_proc(cns, tsne = FALSE)
#' clu_pure(cns, save_dir = save_dir_cr)
#' keep_res <- max_pure(cns, save_dir = save_dir_cr)
#'
#'
max_pure <- function(save_path = paste0(getwd(),
                                        "/outs/all_celltypes/tables/cluster_purity_data"),
                     exponent = 22,
                     pure_thres = 0.99){
  files <- list.files(save_path, pattern = "purity.csv")
  myfiles <- lapply(paste(save_path, files, sep = "/"), read.csv)

  for (j in 1:length(myfiles)) {
    dat <- myfiles[[j]]
    clu_count <- length(levels(as.factor(dat$cluster)))
    clu_pur_mean <- mean(dat$purity)
    pure_measure <- (clu_pur_mean^exponent)/clu_count
    df_temp <- data.frame(clu_res <- files[j],
                          mean_pur = clu_pur_mean,
                          clu_count = clu_count,
                          pure_measure =  pure_measure)


    if(j == 1){
      keep_df <- df_temp
      keep_res <- files[j]
      keep_num_clu<- clu_count
    }else{
      keep_df <- rbind(keep_df, df_temp)
      if(clu_pur_mean >  pure_thres & clu_count > keep_num_clu){
        keep_res <- files[j]
        keep_num_clu<- clu_count
      }
    }
  }

  keep_df$keep_res <- rep(keep_res, nrow(keep_df))
  keep_df$keep_num_clu <- rep(keep_num_clu, nrow(keep_df))


  keep_resolution <- keep_df$keep_res[1]
  keep_resolution <- strsplit(keep_resolution, "_clu")[[1]][1]

  print(paste0("The chosen cluster resolution for downstream analysis is ", keep_resolution))

  #keep_resolution is important for annotation and sub-setting.
  save_summary <- paste0(save_path, "/summary_purity_stats")
  dir.create(save_summary)
  write.csv(keep_df, paste0(save_summary, "/summary_purity_stats.csv"))
  return(keep_resolution)
}







