#' max_pure
#' @description
#' This function uses cluster purity measurements at different resoution and looks
#' for a resolution suitable for a rough celltype annotattion.
#'
#' @param save_dir path to main working directory. Use the same folder as for
#' running cluster_purity.R and the function will find the appropriate files.
#' .csv files. Default is "getwd()/outs/all_celltypes/tables/cluster_purity_data".
#' @param weight_factor this is used to weigh the cluster purity measure. A quotient of
#' cluster purity * weight_factor / cluster count is used to calculate an overall pure
#' measure.
#' @param pure_thres Threshold of cluster purity to be used. The function is looking
#' for the maximum number of clusters over this purity threshold.
#' @param second_thres threshold of the calculated pure measure. Default is 1.
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
#' keep_res <- max_pure(save_dir = save_dir_cr)
#'
#'
max_pure <- function(dir_lab = "all_celltypes",
                     save_dir = getwd(),
                     weight_factor = 22,
                     pure_thres = 0.96,
                     second_thres = 1){
  save_path <- paste0(save_dir,
                      "/outs/",
                      dir_lab,
                      "/tables/cluster_purity_data")
  files <- list.files(save_path, pattern = "purity.csv")
  myfiles <- lapply(paste(save_path, files, sep = "/"), read.csv)

  for (j in 1:length(myfiles)) {
    dat <- myfiles[[j]]
    clu_count <- length(levels(as.factor(dat$cluster)))
    clu_pur_mean <- mean(dat$purity)
    pure_measure <- (clu_pur_mean*weight_factor)/clu_count
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
      if(clu_pur_mean >  pure_thres &
         clu_count > keep_num_clu &
         pure_measure > second_thres) {
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
  if(dir.exists(save_summary) == FALSE){
    dir.create(save_summary)
  }

  write.csv(keep_df, paste0(save_summary, "/summary_purity_stats.csv"))
  return(keep_resolution)
}







