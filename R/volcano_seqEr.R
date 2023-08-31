

volcano_seqEr <- function(dge_res,
                          clu_col = "cluster_id",
                          condition_label = "condition",
                          mode = "overall", # or "clusterwise"
                          save_dir = getwd(),
                          conditions,
                          plotheight = 6,
                          plotwidth = 8){


    if(mode == "overall"){

      cell_levels <- levels(as.factor(dge_res[[clu_col]]))
      dge_res$celltype <- sub("_[^_]+$", "", dge_res[[clu_col]])



        for(j in 1: length(conditions)){
            curr_cond <- cond_levels[j]
            curr_cond_data <- subset(dge_res, dge_res$condition == curr_cond)

            curr_cond_levels <- levels(as.factor(curr_cond_data$cluster))

            if(length(curr_cond_levels) == 2 & min(curr_cond_data$avg_log2FC) > 0){
              proceed = TRUE
              change_logFC <- TRUE
              curr_cond_data$avg_log_2FC_seq <- as.numeric(ifelse(curr_cond_data$cluster == curr_cond_levels[1],
                                                  paste(curr_cond_data$avg_log2FC),
                                                  paste(curr_cond_data$avg_log2FC *(-1))))



            }else if(length(curr_cond_levels) == 1 & min(curr_cond_data$avg_log2FC) < 0) {
              proceed = TRUE
              change_logFC <- FALSE
            }

            if(proceed == TRUE){
                cell_type_lev <- levels(as.factor(curr_cond_data$celltype))
                clu_col_length <- length(cell_type_lev)

                col_df <- data.frame(names = cell_type_lev,
                                     cols = "NA")

                for(k in 1: nrow(col_df)){
                  col_df[k,2] <- colour_palette()[k]
                }


                to_merge_df <- data.frame(names = curr_cond_data$celltype)

                merged_df <- merge(to_merge_df, col_df, by = "names", x.all = TRUE)


                keyvals.col <- merged_df$cols

                names(keyvals.col) <- merged_df$names


                #Plot
                top5 <- curr_cond_data %>% top_n(5, abs(avg_log2FC))
                for (o in 1:length(curr_cond_data$X.1)){
                  ifelse(curr_cond_data$X.1[o] %in% top5$X.1 == TRUE,
                         curr_cond_data$X.1[o] <- curr_cond_data$gene[o],
                         curr_cond_data$X.1[o] <- curr_cond_data$X.1[o])
                }


                p1 <-EnhancedVolcano(curr_cond_data,
                                     lab = curr_cond_data$gene,
                                     x = ifelse(change_logFC == FALSE,
                                                paste0("avg_log2FC"),
                                                names(curr_cond_data)[ncol(curr_cond_data)]) ,
                                     y = "p_val_adj",
                                     labSize = 3.5,
                                     FCcutoff = 0.8,
                                     pCutoff = 0.05,
                                     raster = F,
                                     #shapeCustom = keyvals.shape,
                                     colCustom = keyvals.col,
                                     selectLab = top5$gene,
                                     boxedLabels = T,
                                     drawConnectors = T,
                                     labFace = "bold",
                                     title = curr_cond,
                                     subtitle = paste0("<-- ",
                                                       curr_cond_levels[2],
                                                       " vs. ",
                                                       curr_cond_levels[1],
                                                       " -->"),
                                     legendLabSize = 12,
                                     legendIconSize = 4,
                                     legendPosition = "right")



                save_plot_dir <- paste0(save_dir,
                                        "/outs/summary_figures/conditions/volcano")
                if(dir.exists(save_plot_dir) == FALSE){
                  dir.create(save_plot_dir, recursive = TRUE)
                }

                pdf(paste0(save_plot_dir, "/", curr_cond, "_", mode, ".pdf"),
                    height=plotheight, width = plotwidth)
                print(p1)
                dev.off()

                proceed <-  FALSE
          }



            if(length(curr_cond_levels) > 2){
              print("Consider pairwise comparison if there are more than 3 factor levels
                  in the cluster column or re-run differential gene expression analysis
                  allowing for negative log2FC if only one level is presented in the
                  data frame. Only positive log FC will be considered here.")

              for(p in 1: length(curr_cond_levels)){
                act_clu <- curr_cond_levels[p]
                subs_data <- subset(curr_cond_data, curr_cond_data$cluster == act_clu)

                cell_type_lev <- levels(as.factor(subs_data$celltype))
                clu_col_length <- length(cell_type_lev)

                col_df <- data.frame(names = cell_type_lev,
                                     cols = "NA")

                for(k in 1: nrow(col_df)){
                  col_df[k,2] <- colour_palette()[k]
                }

                to_merge_df <- data.frame(names = subs_data$celltype)

                merged_df <- merge(to_merge_df, col_df, by = "names", x.all = TRUE)


                keyvals.col <- merged_df$cols

                names(keyvals.col) <- merged_df$names


                #Plot
                top5 <- subs_data %>% top_n(5, abs(avg_log2FC))
                for (o in 1:length(subs_data$X.1)){
                  ifelse(subs_data$X.1[o] %in% top5$X.1 == TRUE,
                         subs_data$X.1[o] <- subs_data$gene[o],
                         subs_data$X.1[o] <- subs_data$X.1[o])
                }


                p1 <- EnhancedVolcano(subs_data,
                                     lab = subs_data$gene,
                                     x = "avg_log2FC",
                                     y = "p_val_adj",
                                     labSize = 3.5,
                                     FCcutoff = 0.8,
                                     pCutoff = 0.05,
                                     raster = F,
                                     colCustom = keyvals.col,
                                     selectLab = top5$gene,
                                     boxedLabels = T,
                                     drawConnectors = T,
                                     labFace = "bold",
                                     title = curr_cond,
                                     subtitle = paste0(curr_cond_levels[p],
                                                       " -->"),
                                     legendLabSize = 12,
                                     legendIconSize = 4,
                                     legendPosition = "right")



                save_plot_dir <- paste0(save_dir,
                                        "/outs/summary_figures/conditions/volcano")
                if(dir.exists(save_plot_dir) == FALSE){
                  dir.create(save_plot_dir, recursive = TRUE)
                }

                pdf(paste0(save_plot_dir, "/", curr_cond, "_", curr_cond_levels[p],
                           "_", mode, ".pdf"),
                    height=plotheight, width = plotwidth)
                print(p1)
                dev.off()


              }

            }
          }

        }

  if(mode == "clusterwise"){
    cell_levels <- levels(as.factor(dge_res[[clu_col]]))
    dge_res$celltype <- sub("_[^_]+$", "", dge_res[[clu_col]])

    cell_type_levels <- levels(as.factor(dge_res$celltype))

    for(g in 1: length(cell_type_levels)){
      curr_cell_type <- cell_type_levels[g]
      celltype_data <- subset(dge_res, dge_res$celltype == curr_cell_type)

      for(j in 1: length(conditions)){
        curr_cond <- cond_levels[j]
        curr_cond_data <- subset(celltype_data, celltype_data$condition == curr_cond)

        curr_cond_levels <- levels(as.factor(curr_cond_data$cluster))

        if(length(curr_cond_levels) == 2 & min(curr_cond_data$avg_log2FC) > 0){
          proceed = TRUE
          change_logFC <- TRUE
          curr_cond_data$avg_log_2FC_seq <- as.numeric(ifelse(curr_cond_data$cluster == curr_cond_levels[1],
                                                              paste(curr_cond_data$avg_log2FC),
                                                              paste(curr_cond_data$avg_log2FC *(-1))))


        }else if(length(curr_cond_levels) == 1 & min(curr_cond_data$avg_log2FC) < 0) {
          proceed = TRUE
          change_logFC <- FALSE
        }


        if(proceed == TRUE){
          sub_clu<- levels(as.factor(celltype_data[[clu_col]]))
          sub_clu_length <- length(sub_clu)

          col_df <- data.frame(names = sub_clu,
                               cols = "NA")

          for(k in 1: nrow(col_df)){
            col_df[k,2] <- colour_palette()[k]
          }

          to_merge_df <- data.frame(names = curr_cond_data[[clu_col]])

          merged_df <- merge(to_merge_df, col_df, by = "names", x.all = TRUE)


          keyvals.col <- merged_df$cols

          names(keyvals.col) <- merged_df$names


          #Plot
          top5 <- curr_cond_data %>% top_n(5, abs(avg_log2FC))
          for (o in 1:length(curr_cond_data$X.1)){
            ifelse(curr_cond_data$X.1[o] %in% top5$X.1 == TRUE,
                   curr_cond_data$X.1[o] <- curr_cond_data$gene[o],
                   curr_cond_data$X.1[o] <- curr_cond_data$X.1[o])
          }


          p1 <-EnhancedVolcano(curr_cond_data,
                               lab = curr_cond_data$gene,
                               x = ifelse(change_logFC == FALSE,
                                          paste0("avg_log2FC"),
                                          names(curr_cond_data)[ncol(curr_cond_data)]) ,
                               y = "p_val_adj",
                               labSize = 3.5,
                               FCcutoff = 0.8,
                               pCutoff = 0.05,
                               raster = F,
                               #shapeCustom = keyvals.shape,
                               colCustom = keyvals.col,
                               #selectLab = top5$gene,
                               boxedLabels = T,
                               drawConnectors = T,
                               labFace = "bold",
                               title = curr_cond,
                               subtitle = paste0("<-- ",
                                                 curr_cond_levels[2],
                                                 " vs. ",
                                                 curr_cond_levels[1],
                                                 " -->"),
                               legendLabSize = 12,
                               legendIconSize = 4,
                               legendPosition = "right")



          save_plot_dir <- paste0(save_dir,
                                  "/outs/summary_figures/conditions/volcano")
          if(dir.exists(save_plot_dir) == FALSE){
            dir.create(save_plot_dir, recursive = TRUE)
          }

          pdf(paste0(save_plot_dir, "/", curr_cond, "_", curr_cell_type, "_",
                     mode, ".pdf"),
              height=plotheight, width = plotwidth)
          print(p1)
          dev.off()

          proceed <-  FALSE
        }



        if(length(curr_cond_levels) > 2){
          print("Consider pairwise comparison if there are more than 3 factor levels
                in the cluster column or re-run differential gene expression analysis
                allowing for negative log2FC if only one level is presented in the
                data frame. Only positive log FC will be considered here.")

          for(p in 1: length(curr_cond_levels)){
            act_clu <- curr_cond_levels[p]
            subs_data <- subset(curr_cond_data, curr_cond_data$cluster == act_clu)

            cell_lineage_lev <- levels(as.factor(subs_data[[clu_col]]))
            clu_col_length <- length(cell_lineage_lev)

            col_df <- data.frame(names = cell_lineage_lev,
                                 cols = "NA")

            for(k in 1: nrow(col_df)){
              col_df[k,2] <- colour_palette()[k]
            }

            to_merge_df <- data.frame(names = subs_data[[clu_col]])

            merged_df <- merge(to_merge_df, col_df, by = "names", x.all = TRUE)


            keyvals.col <- merged_df$cols

            names(keyvals.col) <- merged_df$names


            #Plot
            top5 <- subs_data %>% top_n(5, abs(avg_log2FC))
            for (o in 1:length(subs_data$X.1)){
              ifelse(subs_data$X.1[o] %in% top5$X.1 == TRUE,
                     subs_data$X.1[o] <- subs_data$gene[o],
                     subs_data$X.1[o] <- subs_data$X.1[o])
            }


            p1 <- EnhancedVolcano(subs_data,
                                  lab = subs_data$gene,
                                  x = "avg_log2FC",
                                  y = "p_val_adj",
                                  labSize = 3.5,
                                  FCcutoff = 0.8,
                                  pCutoff = 0.05,
                                  raster = F,
                                  #shapeCustom = keyvals.shape,
                                  colCustom = keyvals.col,
                                  selectLab = top5$gene,
                                  boxedLabels = T,
                                  drawConnectors = T,
                                  labFace = "bold",
                                  title = curr_cond,
                                  subtitle = paste0(curr_cond_levels[p],
                                                    " -->"),
                                  legendLabSize = 12,
                                  legendIconSize = 4,
                                  legendPosition = "right")




            save_plot_dir <- paste0(save_dir,
                                    "/outs/summary_figures/conditions/volcano")
            if(dir.exists(save_plot_dir) == FALSE){
              dir.create(save_plot_dir, recursive = TRUE)
            }

            pdf(paste0(save_plot_dir, "/", curr_cond, "_",curr_cell_type, "_",  curr_cond_levels[p],
                       "_", mode, ".pdf"),
                height=plotheight, width = plotwidth)
            print(p1)
            dev.off()
          }
        }
      }
    }
  }
}









