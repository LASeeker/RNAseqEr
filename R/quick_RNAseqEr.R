#' Quick RNAseqEr: Complete workflow for single-cell RNA-seq analysis
#' @description
#' This function runs the complete RNAseqEr workflow from Seurat processing to 
#' final analysis including clustering, annotation, differential expression, 
#' and Shiny app generation.
#'
#' @param seur_obj Seurat object to analyze
#' @param species Species analyzed (default: "Human")
#' @param n_pcs number of principal components for analysis
#' @param res clustering resolutions to test
#' @param select_genes genes to use for scaling
#' @param elbow_dims dimensions for elbow plot
#' @param tsne whether to run tSNE
#' @param col_pattern pattern for clustering column names
#' @param plot_cols colors for plots
#' @param clust_lab whether to show cluster labels
#' @param label_size size of labels
#' @param save_dir directory to save results
#' @param width plot width
#' @param height plot height
#' @param use_reduction dimensional reduction to use
#' @param dir_lab label for output directory
#' @param int_cols columns for differential expression analysis
#' @param only_pos only positive markers
#' @param min_pct minimum percentage of cells
#' @param logfc_threshold log fold change threshold
#' @param fil_pct_1 first percentage filter
#' @param fil_pct_2 second percentage filter
#' @param test_use statistical test to use
#' @param int_cols_pw columns for pairwise analysis
#' @param min_pct_pw minimum percentage for pairwise
#' @param logfc_threshold_pw log fold change threshold for pairwise
#' @param fil_pct_1_pw first percentage filter for pairwise
#' @param fil_pct_2_pw second percentage filter for pairwise
#' @param assay_use assay to use
#' @param ad_pval adjusted p-value threshold
#' @param avg_log average log fold change
#' @param pct_1 first percentage
#' @param pct_2 second percentage
#' @param n_top number of top genes
#' @param use_resol use resolution for heatmaps
#' @param col_names column names for heatmaps
#' @param col_pattern_hm pattern for heatmap columns
#' @param label_hm label heatmaps
#' @param draw_lines draw lines in heatmaps
#' @param min_log2FC_go minimum log2 fold change for GO
#' @param reverse_go reverse for GO analysis
#' @param translate_gene_id_from gene ID format to translate from
#' @param translate_gene_id_to gene ID format to translate to
#' @param org_use organism database
#' @param ontology GO ontology
#' @param pvalue_cutoff_go p-value cutoff for GO
#' @param qvalue_cutoff_go q-value cutoff for GO
#' @param read_able readable gene names
#' @param tissue_ref_annotation tissue reference for annotation
#' @param min_pct_ann minimum percentage for annotation
#' @param logfc_threshold_ann log fold change threshold for annotation
#' @param assay_use_ann assay to use for annotation
#' @param custom_ref_genes_ann use custom reference genes
#' @param custom_gene_list custom gene list
#' @param customclassif custom classification
#' @param plot_reduction reduction for plotting
#' @param plot_label_ann plot labels for annotation
#' @param plot_repel_ann repel labels for annotation
#' @param plot_width_ann plot width for annotation
#' @param plot_height_ann plot height for annotation
#' @param dge_present_ann DGE results present
#' @param dge_ann perform DGE for annotation
#' @param proportion proportion threshold
#' @param cluster_qc_vars variables for cluster QC
#' @param sample_id sample ID column
#' @param cols_milo_design columns for Milo design
#' @param milo_test_fact factors to test with Milo
#' @param sketch whether to use sketched data for large datasets (Seurat 5 feature)
#' @param sketch_ncells number of cells to use for sketching (default: 50000)
#' @param sketch_method sketching method ("LeverageScore" or "Uniform")
#' @param sketched_assay_name name of sketched assay (default: "sketch")
#' @param shiny_name name for the Shiny app (default: "all_celltypes")
#' @param read_file whether to read from file for Shiny app (default: TRUE)
#' @param ext_pattern file extension pattern for Shiny app (default: ".RDS")
#' @param default_1 first default metadata column for Shiny app. If NULL, uses 4th column; if numeric, uses that column index (default: NULL)
#' @param default_2 second default metadata column for Shiny app. If NULL, uses 5th column; if numeric, uses that column index (default: NULL)
#' @param assay_use_shiny assay to use in Shiny app (default: "RNA")
#' @param gex_slot gene expression slots to include in Shiny app (default: c("data", "scale.data", "counts"))
#' @param gene_mapping whether to use gene mapping in Shiny app (default: FALSE)
#' @param default_gene1 first default gene for Shiny app (default: "MALAT1")
#' @param default_gene2 second default gene for Shiny app (default: "GAPDH")
#' @param default_multigene default multi-gene selection for Shiny app (default: NA)
#' @param default_dimred default dimensional reduction for Shiny app (default: c("umap1", "umap2"))
#' @param author author name for Shiny app (default: "My name")
#' @param title publication title for Shiny app (default: "TBC")
#' @param journal journal name for Shiny app (default: "TBC")
#' @param volume journal volume for Shiny app (default: "TBC")
#' @param page journal page for Shiny app (default: "TBC")
#' @param year publication year for Shiny app (default: "2024")
#' @param doi DOI for Shiny app (default: "TBC")
#' @param link link for Shiny app (default: "TBC")
#' @param shiny_title title for Shiny app (default: "My Shiny")
#' @param volcano_mode modes for volcano plots ("overall", "clusterwise", or both) (default: c("overall", "clusterwise"))
#' @param volcano_plotheight height of volcano plots (default: 6)
#' @param volcano_plotwidth width of volcano plots (default: 8)
#'
#' @return Processed Seurat object with all analyses completed
#' @export
#'
#' @examples
#' \dontrun{
#' # Run complete workflow
#' result <- quick_RNAseqEr(cns,
#'                          cluster_qc_vars = c("process_number", "caseNO", "Tissue"),
#'                          sample_id = "uniq_id",
#'                          cols_milo_design = c("uniq_id","Tissue","gender","AgeGroup","caseNO"),
#'                          milo_test_fact = c("Tissue", "gender", "AgeGroup"))
#' 
#' # For large datasets, use sketching
#' result <- quick_RNAseqEr(large_dataset,
#'                          sketch = TRUE,
#'                          sketch_ncells = 100000,
#'                          cluster_qc_vars = c("process_number", "caseNO", "Tissue"),
#'                          sample_id = "uniq_id",
#'                          cols_milo_design = c("uniq_id","Tissue","gender","AgeGroup","caseNO"),
#'                          milo_test_fact = c("Tissue", "gender", "AgeGroup"))
#' }
#'
quick_RNAseqEr <- function(seur_obj,
                           species = "Human",

                           #For Seurat processing
                           n_pcs = 20,
                           nfeatures = 2000,
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

                           #for selecting most appropriate cluster purity
                           weight_factor = 22,
                           pure_thres = 0.96,
                           second_thres = 1,

                           # perform differential gene expression at different
                           #clustering resolutions
                           int_cols = paste0(col_pattern, res),
                           only_pos = TRUE,
                           min_pct = 0.25,
                           logfc_threshold = 0.25,
                           fil_pct_1 = 0.25,
                           fil_pct_2 = 0.6,
                           test_use = "MAST",
                           assay_use = "RNA",


                           #for creating heatmaps
                           use_resol = FALSE,
                           max_diff_threshold = 10,
                           mean_diff_thres = 0.1,

                           #for running gene ontology analysis
                           #gene_list,
                           min_log2FC_go = 0.25,
                           reverse_go = FALSE,
                           translate_gene_id_from = "SYMBOL",
                           translate_gene_id_to = "ENTREZID",
                           org_use = NULL,  # Will be set automatically based on species
                           ontology = "BP",
                           pvalue_cutoff_go = 0.05,
                           qvalue_cutoff_go = 0.05,
                           read_able = TRUE,

                           #specify parameters for annotation using scType
                           tissue_ref_annotation = NULL,  # Will be set automatically based on species and tissue
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
                           plot_height_ann = 6,
                           dge_present_ann = FALSE,
                           dge_ann = TRUE,
                           proportion = 3,

                           #cluster qc
                           cluster_qc_vars,

                           # Milo

                           sample_id,  # use"uniq_id"
                           cols_milo_design, # use c("uniq_id","Tissue","gender","AgeGroup","caseNO")
                           milo_test_fact,  # use c("Tissue", "gender", "AgeGroup")
                           
                           # Sketching parameters for large datasets
                           sketch = FALSE,
                           sketch_ncells = 50000,
                           sketch_method = c("LeverageScore", "Uniform"),
                           sketched_assay_name = "sketch",
                           
                           # Shiny app parameters
                           shiny_name = "all_celltypes",
                           read_file = TRUE,
                           ext_pattern = ".RDS",
                           default_1 = NULL,  # If NULL, will use 4th metadata column; if numeric, will use that column index
                           default_2 = NULL,  # If NULL, will use 5th metadata column; if numeric, will use that column index
                           assay_use_shiny = "RNA",
                           gex_slot = c("data", "scale.data", "counts"),
                           gene_mapping = FALSE,
                           default_gene1 = "MALAT1",
                           default_gene2 = "GAPDH",
                           default_multigene = NA,
                           default_dimred = c("umap1", "umap2"),
                           author = "My name",
                           title = "TBC",
                           journal = "TBC",
                           volume = "TBC",
                           page = "TBC",
                           year = "2024",
                           doi = "TBC",
                           link = "TBC",
                           shiny_title = "My Shiny",
                           
                                  # Volcano plot parameters
       volcano_mode = c("overall", "clusterwise"),
       volcano_plotheight = 6,
       volcano_plotwidth = 8,
       volcano_additional_p_adjust = FALSE,
       volcano_p_adjust_method = "BH",
       
       # Data validation parameters for small datasets
       min_cells_per_cluster = 3,
       min_cells_for_comparison = 5,
       
       # Manuscript generation parameters
       manuscript_generation = FALSE,
       manuscript_title = NULL,
       manuscript_author = NULL,
       manuscript_journal = "bioRxiv",
       manuscript_llm_provider = "openai",
       manuscript_api_key = NULL,
       manuscript_model = "gpt-4",
       manuscript_max_tokens = 4000,
       manuscript_temperature = 0.7,
       manuscript_include_figures = TRUE,
       manuscript_focus_genes = "condition",
       manuscript_min_logfc = 0.5,
       manuscript_max_pvalue = 0.05,
       manuscript_auto_load_env = TRUE
                           
                           #Generating Shiny app
                           ){
  {
    save_dir_save <- save_dir
    
    # Set species-specific parameters
    if(is.null(org_use)) {
      org_use <- get_organism_db(species)
    }
    
    if(is.null(tissue_ref_annotation)) {
      tissue_ref_annotation <- get_tissue_reference(species, "Brain")  # Default to Brain, but could be made configurable
    }
    
    # Validate species support
    if(!is_species_supported(species)) {
      warning(paste("Species", species, "may not be fully supported. Using human defaults as fallback."))
    }
    
    # Analyze dataset size and provide recommendations
    print("Step 0/8: Analyzing dataset size...")
    size_analysis <- analyze_dataset_size(seur_obj, verbose = TRUE)
    
    # Adjust parameters based on dataset size if needed
    if(size_analysis$size_category == "very_small") {
      cat("⚠️  Very small dataset detected. Adjusting parameters for reliability...\n")
      min_cells_per_cluster <- max(min_cells_per_cluster, 3)
      min_cells_for_comparison <- max(min_cells_for_comparison, 5)
    }
  
  
    # Perform standard Seurat processing
  
    print("Step 1/9: Performing Seurat processing...")
    seur_obj <- seurat_proc(seur_obj,
                            n_pcs = n_pcs,
                            res = res,
                            select_genes = select_genes,
                            elbow_dims = elbow_dims,
                            tsne = tsne,
                            sketch = sketch,
                            sketch_ncells = sketch_ncells,
                            sketch_method = sketch_method,
                            sketched_assay_name = sketched_assay_name)
  
    # plot dimensionally reduced data at different clustering resolutions and
    # save to file
  
  
    print("Step 2/9: Plotting dimensional reductions...")
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
  
  
    print("Step 3/9: Calculating cluster purity...")
    clu_pure(seur_obj,
             reduction = use_reduction, #reduction_sil,
             col_pattern = col_pattern,
             plot_cols = plot_cols,
             clust_lab = clust_lab,
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
                               ident_use = keep_res,
                               tissue_ref = tissue_ref_annotation,
                               test_use = test_use,
                               min_pct = min_pct,
                               logfc_threshold = logfc_threshold_ann,
                               assay_use = if(sketch) sketched_assay_name else assay_use,
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
                      test_use = test_use,
                      assay_use = if(sketch) sketched_assay_name else assay_use)
  
    cond_genes <- gen_mark_list(file_dir = paste0(save_dir,
                                                  "/outs/",
                                                  dir_lab,
                                                  "/tables/condition_mark/",
                                                  "RNAseqEr_annotation",
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
    
    # Generate volcano plots for condition markers
    print("Step 6/9: Generating volcano plots for condition markers...")
    
    # Read in the condition marker results for volcano plotting
    cond_mark_dir <- paste0(save_dir, "/outs/", dir_lab, "/tables/condition_mark/RNAseqEr_annotation")
    
    if(dir.exists(cond_mark_dir)) {
      # Generate volcano plots based on specified modes
      for(mode in volcano_mode) {
                 volcano_seqEr(dge_res = cond_genes,
                       clu_col = "cluster",
                       condition_label = "condition",
                       mode = mode,
                       save_dir = save_dir,
                       conditions = milo_test_fact,
                       plotheight = volcano_plotheight,
                       plotwidth = volcano_plotwidth,
                       additional_p_adjust = volcano_additional_p_adjust,
                       p_adjust_method = volcano_p_adjust_method)
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
                                       dir_lab = cell_lineage_name,
                                       assay_use = if(sketch) sketched_assay_name else assay_use,
                                       min_cells_per_cluster = min_cells_per_cluster,
                                       min_cells_for_comparison = min_cells_for_comparison)
  
      pw_mark <- pairwise_dge(seur_obj = curr_srt,
                              int_cols = keep_res_ul,
                              save_dir = save_dir,
                              dir_lab = cell_lineage_name,
                              assay_use = if(sketch) sketched_assay_name else assay_use,
                              min_cells_per_cluster = min_cells_per_cluster,
                              min_cells_for_comparison = min_cells_for_comparison)
  
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
  
            # Save the final clustering plot
      final_cluster_plot <- DimPlot(curr_srt, group.by = "cluster_id", label = TRUE, cols = colour_palette())
      
      # Save the plot
      cluster_plot_dir <- paste0(save_dir, "/outs/", cell_lineage_name, "/plots/final_clustering/")
      if(dir.exists(cluster_plot_dir) == FALSE){
        dir.create(cluster_plot_dir, recursive = TRUE)
      }
      
      pdf(paste0(cluster_plot_dir, "final_clustering_", cell_lineage_name, ".pdf"), 
          width = 10, height = 8)
      print(final_cluster_plot)
      dev.off()
  
  
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
                                  dir_lab = cell_lineage_name,
                                  assay_use = if(sketch) sketched_assay_name else assay_use)
  
  
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
                                      exclude_out_dir = c("data","supplementary_tables"),
                                      comp_mode = comp_modes[v],
                                      conditions = milo_test_fact,
                                      cluster_label = "cell_lineage")
      
             # Generate volcano plots for condensed marklists
       if(nrow(curr_mark) > 0) {
         volcano_seqEr(dge_res = curr_mark,
                       clu_col = "cluster_id",
                       condition_label = "condition",
                       mode = comp_modes[v],
                       save_dir = save_dir,
                       conditions = milo_test_fact,
                       plotheight = volcano_plotheight,
                       plotwidth = volcano_plotwidth,
                       additional_p_adjust = volcano_additional_p_adjust,
                       p_adjust_method = volcano_p_adjust_method)
       }
  
    }
  
  
    #save annotated dataset
    saveRDS(seur_obj, paste0(proc_dat_dir_all, "/", dir_lab, ".RDS"))
  
  
        # build shiny
    
    # Handle default_1 and default_2 intelligently
    if(is.null(default_1)) {
      # Use 4th metadata column if available
      meta_cols <- colnames(seur_obj@meta.data)
      if(length(meta_cols) >= 4) {
        default_1 <- meta_cols[4]
      } else {
        default_1 <- meta_cols[1]  # fallback to first column
      }

    create_shiny(seur_obj,
                 shiny_name = shiny_name,
                 read_file = read_file,
                 file_dir_1 = proc_dat_dir_all,
                 file_dir_2 = proc_dat_dir,
                 ext_pattern = ext_pattern,
                 default_1 = default_1,
                 default_2 = default_2,
                 assay_use = assay_use_shiny,
                 gex_slot = gex_slot,
                 gene_mapping = gene_mapping,
                 default_gene1 = default_gene1,
                 default_gene2 = default_gene2,
                 save_dir = save_dir,
                 default_multigene = default_multigene,
                 default_dimred = default_dimred,
                 author = author,
                 title = title,
                 journal = journal,
                 volume = volume,
                 page = page,
                 year = year,
                 doi = doi,
                 link = link,
                 shiny_title = shiny_title
    )
  
  
    # Generate comprehensive analysis report
    print("Step 8/9: Generating analysis reports...")
    
    # Collect workflow parameters for the report
    workflow_params <- list(
      n_pcs = n_pcs,
      res = res,
      test_use = test_use,
      tissue_ref = tissue_ref_annotation,
      cluster_qc_vars = cluster_qc_vars,
      logfc_threshold = logfc_threshold,
      min_pct = min_pct,
      pvalue_cutoff_go = pvalue_cutoff_go,
      qvalue_cutoff_go = qvalue_cutoff_go,
      keep_res = keep_res,
      pure_thres = 0.96,  # Default from max_pure function
      weight_factor = 22,  # Default from max_pure function
      ad_pval = ad_pval,
      avg_log = avg_log,
      # Volcano plot parameters
      volcano_mode = volcano_mode,
      volcano_additional_p_adjust = volcano_additional_p_adjust,
      volcano_p_adjust_method = volcano_p_adjust_method
    )
    
    # Generate reports
    report_dir <- generate_report(seur_obj = seur_obj,
                                save_dir = save_dir,
                                dir_lab = dir_lab,
                                workflow_params = workflow_params,
                                analysis_date = Sys.Date(),
                                researcher_name = "Researcher",
                                project_title = "RNAseqEr Single-cell Analysis",
                                tissue_type = tissue_ref_annotation,
                                species = "Human")
    
    # Generate manuscript draft if requested
    if(!is.null(manuscript_generation) && manuscript_generation) {
      print("Step 9/9: Generating manuscript draft...")
      
      manuscript_dir <- generate_manuscript(save_dir = save_dir,
                                          dir_lab = dir_lab,
                                          conditions_of_interest = milo_test_fact,
                                          tissue_type = tissue_ref_annotation,
                                          species = species,
                                          project_title = manuscript_title,
                                          researcher_name = manuscript_author,
                                          journal_target = manuscript_journal,
                                          llm_provider = manuscript_llm_provider,
                                          api_key = manuscript_api_key,
                                          model_name = manuscript_model,
                                          max_tokens = manuscript_max_tokens,
                                          temperature = manuscript_temperature,
                                          include_figures = manuscript_include_figures,
                                          focus_on_genes = manuscript_focus_genes,
                                          min_logfc = manuscript_min_logfc,
                                          max_pvalue = manuscript_max_pvalue,
                                          auto_load_env = manuscript_auto_load_env)
    }
    
    # At the end, print a summary
    cat("\n=== RNAseqEr Analysis Complete ===\n")
    cat("Output saved to:", save_dir, "\n")
    cat("Shiny app created in: shiny_app/\n")
    cat("Processed data saved to: outs/data/processed/\n")
    cat("Analysis reports saved to:", report_dir, "\n")
    if(!is.null(manuscript_generation) && manuscript_generation) {
      cat("Manuscript draft saved to:", manuscript_dir, "\n")
    }
    cat("\nReports generated:")
    cat("- methods_report.md (Methods section for papers)")
    cat("- decision_justifications.md (Justifications for analysis decisions)")
    cat("- complete_analysis_report.md (Combined report)")
    if(!is.null(manuscript_generation) && manuscript_generation) {
      cat("- manuscript_draft.md (LLM-generated manuscript draft)")
      cat("- figure_selection.md (Selected plots for figures)")
      cat("- summary_statistics.md (Key statistics and results)")
    }
  
    return(seur_obj)
  }
