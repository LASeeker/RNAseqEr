
quick_RNAseqEr <- function(seur_obj,

                           #For Seurat processing
                           n_pcs = 20,
                           res = c(0.005, 0.01, 0.04, 0.05,
                                   seq(from = 0.1, to = 1.5, by = 0.1)),
                           select_genes = rownames(seur_obj),
                           elbow_dims = 50,
                           tsne = FALSE,

                           #for plotting dimensional reduction
                           col_pattern = "RNA_snn_res.",
                           plot_cols = colour_palette(),
                           clust_lab = TRUE,
                           label_size = 8,
                           save_dir = getwd(),
                           width=7,
                           height=5,
                           use_reduction = "umap",
                           dir_lab = "all_celltypes",

                           # perform differential gene expression at different
                           #clustering resolutions
                           int_cols = paste0(col_pattern, res),
                           only_pos = TRUE,
                           min_pct = 0.25,
                           logfc_threshold = 0.25,
                           fil_pct_1 = 0.25,
                           fil_pct_2 = 0.6,
                           test_use = "MAST",

                           # pairwise comparison of clusters at one or more
                           # resolutions of interest
                           int_cols_pw = "RNA_snn_res.1",
                           min_pct_pw = 0.25,
                           logfc_threshold_pw = 0.25,
                           fil_pct_1_pw = 0.25,
                           fil_pct_2_pw = 0.1,
                           assay_use = "RNA",

                           # compile list of interesting genes based on top x
                           # hits per cluster/group of interest using DGE
                           # results (pairwise and not)

                           ad_pval = 0.05,
                           avg_log = 1.2,
                           pct_1 = 0.25,
                           pct_2 = 0.6,
                           n_top = 10,

                           #for creating heatmaps
                           use_resol = TRUE,
                           col_names = NA,
                           col_pattern_hm = "RNA_snn_res",
                           label_hm = TRUE,
                           draw_lines = FALSE,

                           #for running gene ontology analysis
                           #gene_list,
                           min_log2FC_go = 0.25,
                           reverse_go = FALSE,
                           translate_gene_id_from = "SYMBOL",
                           translate_gene_id_to = "ENTREZID",
                           org_use = "org.Hs.eg.db",
                           ontology = "BP",
                           pvalue_cutoff_go = 0.05,
                           qvalue_cutoff_go = 0.05,
                           read_able = TRUE,

                           #specify parameters for annotation using scType
                           tissue_ref_annotation = "Brain",
                           min_pct_ann = 0.25,
                           logfc_threshold_ann = 0.8,
                           assay_use_ann = "RNA",
                           custom_ref_genes_ann = FALSE,
                           custom_gene_list = NA,
                           customclassif = "ScType_annotation",
                           plot_reduction = "UMAP",
                           plot_label_ann = TRUE,
                           plot_repel_ann = TRUE,
                           plot_width_ann = 8,
                           plot_height_ann = 8,
                           dge_present_ann = FALSE,
                           dge_ann = TRUE,
                           proportion = 3,

                           #cluster qc
                           cluster_qc_vars, # use c("process_number", "caseNO", "Tissue"),

                           # Milo
                           sample_id,  # use"uniq_id"
                           cols_milo_design, # use c("uniq_id","Tissue","gender","AgeGroup","caseNO")
                           milo_test_fact  # use c("Tissue", "gender", "AgeGroup")
                           
                           #Generating Shiny app
                           
                           
                           
                           ){
  #tryCatch({
    save_dir_save <- save_dir
  
  
  
    # Perform standard Seurat processing
  
    seur_obj <- seurat_proc(seur_obj,
                            n_pcs = n_pcs,
                            res = res,
                            select_genes = select_genes,
                            elbow_dims = elbow_dims,
                            tsne = tsne)
  
    # plot dimensionally reduced data at different clustering resolutions and
    # save to file
  
  
    plot_list(seur_obj = seur_obj,
              col_pattern = col_pattern,
              plot_cols = plot_cols ,
              clust_lab = clust_lab,
              label_size = label_size,
              save_dir = save_dir,
              width=width,
              height=height,
              use_reduction = use_reduction)
  
  
    # Calculate and plot cluster purity measures
  
  
    clu_pure(seur_obj,
             reduction = use_reduction, #reduction_sil,
             col_pattern = col_pattern,
             plot_cols = plot_cols,
             clust_lab = clust_lab_sil,
             label_size = label_size,
             save_dir = save_dir, #clu_pur_dir,
             width=7,
             height=5)
  
  
    #read in purity measures to see which cluster resolution is the one with the
    #largest number of clusters and largest purity measure
  
    keep_res <- max_pure(save_dir = save_dir)
  
    Idents(seur_obj) <- keep_res
  
  
    ## Annotation
    seur_obj <- annotate_seqEr(seur_obj,
                               dir_lab = dir_lab,
                               ident_use = keep_resolution,
                               tissue_ref = tissue_ref_annotation,
                               test_use = test_use,
                               min_pct = min_pct,
                               logfc_threshold = logfc_threshold_ann,
                               assay_use = assay_use,
                               save_dir = save_dir,
                               custom_ref_genes = custom_ref_genes_ann,
                               custom_gene_list,
                               customclassif = customclassif,
                               plot_reduction = use_reduction,
                               plot_label = plot_label_ann,
                               plot_repel = plot_repel_ann,
                               plot_width = plot_width_ann,
                               plot_height = plot_height_ann,
                               dge_present = dge_present_ann,
                               dge = dge_ann,
                               proportion = proportion,
                               colours = plot_cols)
  
  
    if(save_dir != save_dir_save){
      save_dir <- save_dir_save
    }
  
  
    #run differential gene expression for chosen clustering resolution
    Idents(seur_obj) <- "RNAseqEr_annotation"
  
    dge_dir <- paste0(save_dir, "/outs/", dir_lab, "/tables/DGE/broad_celltype_markers")
  
    clu_mark <- gen_mark_list(file_dir = dge_dir)
  
    int_genes <- unique(clu_mark$gene)
  
    hm_dir <- paste0(save_dir,
                     "/outs/",
                    dir_lab,
                   "/plots/",
                   "/heatmaps/RNAseqEr_annot")
  
  
    dir.create(hm_dir, recursive = TRUE)
  
  
    df_seqEr <- heatmap_seqEr(seur_obj,
                              use_resol = FALSE,
                              col_names = "RNAseqEr_annotation",
                              int_genes = int_genes,
                              save_dir = hm_dir)
  
  
    # Do cluster QC
  
    seur_obj <- cluster_qc(seur_obj = seur_obj,
                       cluster_col = "RNAseqEr_annotation",
                       vars = cluster_qc_vars,
                       dir_lab = dir_lab)
  
    DimPlot(seur_obj, group.by = "RNAseqEr_annotation", label = TRUE) +NoLegend()
  
    seur_obj <- remove_clu(seur_obj = seur_obj,
                           cluster_col = "RNAseqEr_annotation",
                           save_dir = save_dir,
                           dir_lab = dir_lab)
  
    DimPlot(seur_obj, group.by = "RNAseqEr_annotation", label = TRUE) +NoLegend()
  
  
    # Test if abundancy differs with variables of technical or biological interest
    milo_obj <- abundance_test(seur_obj = seur_obj,
                               cluster_col = "RNAseqEr_annotation",
                               sample_id,
                               cols_interest = cols_milo_design,
                               test_factors = milo_test_fact,
                               dir_lab = dir_lab)
  
  
  
  
    find_cond_markers(seur_obj,
                      int_cols = milo_test_fact,
                      cluster_id = "RNAseqEr_annotation",
                      dir_lab = dir_lab,
                      only_pos = only_pos,
                      min_pct = min_pct,
                      logfc_threshold = logfc_threshold,
                      fil_pct_1 = fil_pct_1,
                      fil_pct_2 = fil_pct_2,
                      save_dir = save_dir_save,
                      test_use = test_use)
  
    cond_genes <- gen_mark_list(file_dir = paste0(save_dir,
                                                  "/outs/",
                                                  dir_lab,
                                                  "/tables/condition_mark/",
                                                  cluster_id,
                                                  "/clusterwise"),
                                condition = TRUE,
                                test_cond = milo_test_fact)
  
    for(g in 1: length(milo_test_fact)){
      curr_fact_name <- milo_test_fact[g]
      if(curr_fact_name %in% cond_genes$condition){
        sub_dat <- subset(cond_genes, cond_genes$condition == curr_fact_name)
        sub_dat <- sub_dat[!duplicated(sub_dat$gene),]
  
        # Save vln plots to file
        save_vln(seur_obj,
               marker_list = sub_dat$gene,
               save_label = paste0(milo_test_fact[g], "_mark"),
               condition_label = "Condition",
               split_by = milo_test_fact[g],
               group_by = "RNAseqEr_annotation",
               plotheight = 50,
               save_dir = save_dir)
  
        save_feat_plots(seur_obj,
                        sub_dat$gene,
                        dir_lab = "all_celltypes",
                        save_label = paste0(milo_test_fact[g], "_mark"),
                        condition_label = "Condition",
                        save_dir = getwd(),
                        numb_genes = 3,
                        plotheight = 18,
                        plotwidth = 18,
                        split_by = milo_test_fact[g],
                        n_col = 4)
      }
    }
  
    # Now that we have a resolution for the dataset that enabled a preliminary annotatoion,
    # let's subset for the main clusters (as in cell lineages) and look for
    # finer clusters that make sense biologically
  
  
    subset_RNAseqEr(seur_obj,
                    subset_column = "RNAseqEr_annotation",
                    save_dir = save_dir)
  
  
  
    # process cell lineage datasets
  
    save_path <- paste0(save_dir, "/outs/data")
    cl_dat_name <- list.files(save_path)
  
    for(y in 3: length(cl_dat_name)){
      cur_name <- cl_dat_name[y]
      curr_srt <- readRDS(paste0(save_path, "/", cur_name))
  
      cell_lineage_name <- strsplit(cur_name, ".R")[[1]][1]
  
      curr_srt <- seurat_proc(curr_srt,
                             n_pcs = n_pcs,
                             res = res,
                             elbow_dims = elbow_dims,
                             tsne = tsne,
                             dir_lab = cell_lineage_name)
  
      plot_list(seur_obj = curr_srt,
                col_pattern = col_pattern,
                plot_cols = plot_cols ,
                clust_lab = clust_lab,
                label_size = label_size,
                save_dir = save_dir,
                width=width,
                height=height,
                use_reduction = use_reduction,
                dir_lab = cell_lineage_name)
  
      keep_res <- select_res(curr_srt)
      keep_res_ul <- unlist(keep_res)
  
  
  
      all_res_mark <- int_res_all_mark(seur_obj = curr_srt,
                                       int_cols = keep_res_ul,
                                       save_dir = save_dir,
                                       dir_lab = cell_lineage_name)
  
      pw_mark <- pairwise_dge(seur_obj = curr_srt,
                              int_cols = keep_res_ul,
                              save_dir = save_dir,
                              dir_lab = cell_lineage_name)
  
      save_path_o <- paste0(save_dir,
                          "/outs/",
                          cell_lineage_name,
                          "/tables/cluster_marker/overall"
      )
  
      save_path_pw <- paste0(save_dir,
                            "/outs/",
                            cell_lineage_name,
                            "/tables/cluster_marker/pairwise"
      )
  
  
      clu_mark <- gen_mark_list(file_dir = save_path_o)
  
      clu_mark_pw <- gen_mark_list(file_dir = save_path_pw,
                                   pairwise = TRUE)
  
      int_genes <- c(clu_mark$gene, clu_mark_pw$gene)
      int_genes <- unique(int_genes)
  
  
      hm_dir <- paste0(save_dir,
                       "/outs/",
                       cell_lineage_name,
                     "/plots/heatmaps")
  
  
      if(dir.exists(hm_dir) == FALSE){
        dir.create(hm_dir, recursive = TRUE)
      }
  
  
      df_seqEr <- heatmap_seqEr(curr_srt,
                                use_resol = FALSE,
                                col_names = keep_res_ul,
                                int_genes = int_genes,
                                save_dir = hm_dir,
                                max_diff_threshold = 10,
                                mean_diff_thres = 0.1)
  
      df_seqEr$res <- as.numeric(sapply(strsplit(df_seqEr$resolution, col_pattern), "[", 2))
  
      sub_df <- subset(df_seqEr, df_seqEr$cluster_sim_score == min(df_seqEr$cluster_sim_score))
      sub_df <- subset(sub_df, sub_df$res == max(sub_df$res))
  
      keep_res <- sub_df$resolution
  
  
      Idents(curr_srt) <- keep_res
  
      curr_srt$cluster_id <- paste0(cell_lineage_name, "_", curr_srt@meta.data[[keep_res]])
  
      DimPlot(curr_srt, group.by = "cluster_id", label = TRUE, cols = colour_palette())
  
      #TO DO SAVE THE PLOT ABOVE
  
  
      curr_srt <- cluster_qc(seur_obj = curr_srt,
                         cluster_col = "cluster_id",
                         vars = cluster_qc_vars,
                         dir_lab = cell_lineage_name)
  
  
      milo_obj <- abundance_test(seur_obj = curr_srt,
                                 cluster_col = "cluster_id",
                                 sample_id = sample_id,
                                 cols_interest = cols_milo_design,
                                 test_factors = milo_test_fact,
                                 dir_lab = cell_lineage_name)
  
  
  
      RNAseqEr::find_cond_markers(curr_srt,
                                  int_cols = milo_test_fact,
                                  cluster_id = "cluster_id",
                                  dir_lab = cell_lineage_name)
  
  
      cond_genes <- gen_mark_list(file_dir = paste0(save_dir,
                                                    "/outs/",
                                                    cell_lineage_name,
                                                  "/tables/",
                                                  "condition_mark/" ,
                                                  "cluster_id/" ,
                                                  "clusterwise"),
                                  condition = TRUE,
                                  test_cond = milo_test_fact)
  
  
      for(q in 1: length(milo_test_fact)){
        mark <- subset(cond_genes,
                       cond_genes$condition == milo_test_fact[q])
        if(nrow(mark) > 0){
          mark_red <- mark[!duplicated(mark$gene),]
  
          save_vln(curr_srt,
                   dir_lab = cell_lineage_name,
                   marker_list = mark_red$gene,
                   save_label = paste0(milo_test_fact[q], "_mark_cw"),
                   split_by = milo_test_fact[q],
                   group_by = "cluster_id",
                   plotheight = 50,
                   condition_label = "Condition")
  
          save_feat_plots(curr_srt,
                          marker_list = mark_red$gene,
                          dir_lab = cell_lineage_name,
                          save_label = paste0(milo_test_fact[q], "_mark_cw"),
                          split_by = milo_test_fact[q],
                          plotheight = 30,
                          condition_label = "Condition")
  
  
        }
      }
  
  
      # Overall
      cond_genes <- gen_mark_list(file_dir = paste0(save_dir,
                                                    "/outs/",
                                                    cell_lineage_name,
                                                    "/tables/",
                                                    "condition_mark/" ,
                                                    "cluster_id/" ,
                                                    "overall"),
                                  condition = TRUE,
                                  test_cond = milo_test_fact)
  
  
      for(q in 1: length(milo_test_fact)){
        mark <- subset(cond_genes,
                       cond_genes$condition == milo_test_fact[q])
        if(nrow(mark) > 0){
          mark_red <- mark[!duplicated(mark$gene),]
  
          save_vln(curr_srt,
                   condition_label = "Condition",
                   marker_list = mark_red$gene,
                   dir_lab = cell_lineage_name,
                   save_label = paste0(milo_test_fact[q], "_mark_ov"),
                   split_by = milo_test_fact[q],
                   group_by = "cluster_id",
                   plotheight = 50)
  
          save_feat_plots(curr_srt,
                          condition_label = "Condition",
                          marker_list = mark_red$gene,
                          dir_lab = cell_lineage_name,
                          save_label = paste0(milo_test_fact[q], "_mark_ov"),
                          split_by = milo_test_fact[q],
                          plotheight = 30)
  
  
          # Run GO
          clu_levels <- levels(as.factor(mark$cluster))
          for(u in 1: length(clu_levels)){
            subs_dat <- subset(mark, mark$cluster == clu_levels[u])
            curr_level <- clu_levels[u]
            
            go_results <- perform_go(curr_srt,
                                     gene_list = mark,
                                     min_log2FC = 0.25,
                                     reverse = FALSE,
                                     translate_gene_id_from = "SYMBOL",
                                     translate_gene_id_to = "ENTREZID",
                                     org_use = "org.Hs.eg.db",
                                     ontology = "BP",
                                     pvalue_cutoff = 0.05,
                                     qvalue_cutoff = 0.05,
                                     read_able = TRUE)
  
            go_plot <- dotplot(go_results)
  
            #go_save_dir <- paste0(save_dir, "/outs/", cell_lineage_name,
            #                      "/plots/gene_ontology/condition")
  
            #if(dir.exists(go_save_dir) == FALSE){
            #  dir.create(go_save_dir, recursive = TRUE)
            #}
  
  
            #pdf(paste0(go_save_dir, "/", curr_level, "_go.pdf"),
            #    height=plotheight, width = plotwidth)
            #print(go_plot)
            #dev.off()
  
            #print(go_plot)
  
          # THe below needs to be moved up into the loop
  
          proc_dat_dir <- paste0(save_dir, "/outs/data/processed")
  
          if(dir.exists(proc_dat_dir) == FALSE){
            dir.create(proc_dat_dir)
          }
          saveRDS(curr_srt, paste0(proc_dat_dir, "/", cell_lineage_name, ".RDS"))
  
  
          }
  
  
        }
    }
  
  
  
    seur_obj <- fine_annotate(seur_obj,
                              file_dir = proc_dat_dir)
  
    proc_dat_dir_all <- paste0(save_dir, "/outs/data/processed/all_data")
    if(dir.exists(proc_dat_dir_all) == FALSE){
      dir.create(proc_dat_dir_all)
    }
  
  
    #MEGA DOTPLOT
  
  
  
    ## COLLECT INFO ON CONDITION
  
    comp_modes <- c("overall", "clusterwise")
  
    for(v in 1:length(comp_modes)){
      curr_mark <- condense_marklists(save_dir = save_dir,
                                      exclude_out_dir = c("data","supplemenraty_tables"),
                                      comp_mode = comp_modes[v],
                                      conditions = milo_test_fact,
                                      cluster_label = "cell_lineage")
      ## ADD here that curr_mark list is plotted as volcano
  
    }
  
  
  
  
  
  
    #save annotated dataset
    saveRDS(seur_obj, paste0(proc_dat_dir_all, "/", dir_lab, ".RDS"))
  
  
    # build shiny
  
  
    create_shiny(seur_obj,
                 shiny_name = "all_celltypes",
                 read_file = TRUE,
                 file_dir_1 = proc_dat_dir_all,
                 file_dir_2 = proc_dat_dir,
                 ext_pattern = ".RDS",
                 default_1 = "AgeGroup",
                 default_2 = "Tissue",
                 assay_use = "RNA",
                 gex_slot = c("data", "scale.data", "counts"),
                 gene_mapping = FALSE,
                 default_gene1 = "MALAT1",
                 default_gene2 = "GAPDH",
                 save_dir = getwd(),
                 default_multigene = NA,
                 default_dimred = c("umap1", "umap2"),
                 author = "TBC",
                 title = "TBC",
                 journal = "TBC",
                 volume  = "TBC",
                 page    = "TBC",
                 year    = "TBC",
                 doi     = "TBC",
                 link    = "TBC",
                 shiny_title = "My Shiny"
    )
  
  
    return(seur_obj)
    } 
    
    #error = function(e) {
    #  message("Error in quick_RNAseqEr: ", e$message)
     # # Save partial results if possible
    #  return(NULL)
    #  }, warning = function(w) {
    #    message("Warning in quick_RNAseqEr: ", w$message)
    #  })
  }
