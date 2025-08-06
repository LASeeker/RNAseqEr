#' Fast RNAseqER analysis of pre-processed single cell/ nuclei datasets
#' @description
#' This function takes a quality controlled dataset (gene and cell QC) with optional
#' other corrections (ambient RNA removel, batch correction, integration) performed
#' and produces results including documentation, seurat opbjects for further analysis
#' and a shiny app.
#'
#'
#' @param seur_obj quality controlled Seurat object
#' @param n_pcs number of Principal components used for dimensional reduction
#' and clustering
#' @param res resolutions that are being used for clustering. The default are the
#' following: 0.005, 0.010, 0.040, 0.050, 0.100, 0.200, 0.300, 0.400, 0.500,
#' 0.600, 0.700, 0.800, 0.900, 1.000, 1.100, 1.200, 1.300, 1.400, 1.500
#' @param select_genes which genes should be used for scaling. The default is
#' all (rownames(seur_obj))
#' @param elbow_dims number of dimensions plotted in Elbow plot to select principal
#' components for dimensional reductions and clustering. The default is 50.
#' @param tsne TRUE/FALSE whether or not TSNE should be performed. default is
#' FALSE.
#' @param col_pattern column pattern used to detect meta data columns that contain
#' cluster information. The default is "RNA_snn_res.".
#' @param plot_cols string of colours used in plots. The default is RNAseqEr's
#' colour_palette()
#' @param clust_lab TRUE/FALSE whether or not cluster labels should be displayed in plot.
#' The default is TRUE.
#' @param label_size size of cluster labels (default 8)
#' @param save_dir path to root directory for analysis.
#' @param width plot width (7)
#' @param height plot height (5)
#' @param use_reduction which dimensiona reduction should be used. Default is "umap".
#' @param dir_lab label for complete dataset. Default is "all_celltype". This can
#' be freely chosen and will be used to generate file paths.
#' @param int_cols columns of interest for the analysis. This could for example
#' be the metadata column names corresponding to treatment vs. control identity,
#' age group, tissue region, KO vs. control etc.
#' @param only_pos TRUE/FALSE whether only positive values should be considered
#' in differential gene expression analyses. The default is TRUE.
#' @param min_pct Minimum proportion of cells within the cluster/group of interest expressing
#' the gene for gene to appear in the differential gene. Default is 0.25.
#' @param logfc_threshold Minimum log2PC for gene to appear in the differential gene
#' expression analysis result. The default is 0.25.
#' @param fil_pct_1 Filter parameter for the proportion of cells WITHIN the
#' cluster/condition of interest expressing a certain genes of the differential gene
#' expression analysis result. Default is 0.25 and thus does not filter.
#' @param fil_pct_2 Filter parameter for the proportion of cells OUTSIDE the
#' cluster/condition of interest expressing a certain genes of the differential gene
#' expression analysis result. Default is 0.6.
#' @param test_use test to use for differential gene expression analysis. The
#' default is "MAST" but all other methods implemented within Seurat work
#' as alternatives.
#' @param assay_use assay to use for differential gene expression analysis.
#' The default is "RNA"
#' @param use_resol TRUE/FALSE, Default is FALS
#' @param min_log2FC_go min log2 fold change that a gene needs to pass to be
#' considered in a gene ontology analysis. Default is 0.25
#' @param reverse_go TRUE/FALSE EXPLAIN. Default is FALSE.
#' @param translate_gene_id_from from which gene annotation should be translated.
#' The default is "SYMBOL".
#' @param translate_gene_id_to To which gene annotation should be translated. The
#' default is "ENTREZID"
#' @param org_use Organism used. The default is human ("org.Hs.eg.db")
#' @param ontology What kind of gene ontology analysis is to be performed. The
#' default is "BP" for biological process, but "MF" for molecular function is
#' an option as well.
#' @param pvalue_cutoff_go max p-value for a GO result to be considered statistically
#' significant.
#' @param qvalue_cutoff_gomax max q-value for a GO result to be considered statistically
#' significant.
#' @param read_able TRUE/FALSE whether gene names should be readable. The default is
#' TRUE.
#' @param tissue_ref_annotation Reference tissue for annotation. The default is
#' "Brain" and all other scType reference tissues work.
#' @param min_pct_ann Minimum pct1 for differential gene expression analysis that will
#' be used for the cell type annotation.
#' @param logfc_threshold_ann Minimum log fold change threshold for differential
#' gene expression analysis that will be used for the cell type annotation.
#' @param assay_use_ann assay used for cell type annotation. The default is "RNA"
#' @param custom_ref_genes_ann TRUE/FALSE whether list of genes for the annotation
#' is provided. Default is FALSE
#' @param custom_gene_list If custom_ref_genes_ann == TRUEa reference gene list
#' has to be provided. Default is that the opposite and value is set to NA
#' @param customclassif name for classification used. Default is "ScType_annotation"
#' @param plot_reduction reduction used for plotting annotation results.
#' Default is "plot_reduction"
#' @param plot_label_ann TRUE/FALSE whether labels should be plotted. Default is
#' TRUE.
#' @param plot_repel_ann TRUE/FALSE whether labels should be repelled Default is
#' TRUE.
#' @param plot_width_ann plot width of annotated dataset. Default = 8
#' @param plot_height_ann plot height of annotated dataset. Default = 6.
#' @param dge_present_ann TRUE/FALSE whether differential gene expression results
#' are already present for example because the funcltion has been run previously.
#' It is important that it has been run using the same settings up to the point of
#' annotation. If so, setting the default FALSE value to TRUE can speed up the
#' annotation process.
#' @param dge_ann TRUE/FALSE whther differential gene expression should be run.
#' Default is TRUE.
#' @param proportion threshold proportion for annotation. Default is set to 3 which
#' means that the number differentially expressed genes in a cluster must be
#' larger than a third of the total number of genes that define a cell type in scTypoe to
#' be considered that cell type. Increasing the number makes annotation less stringent,
#' decreasing it makes it more stringent.
#' @param cluster_qc_vars String of variables to test for cluster quality control.
#' for example c("process_number", "caseNO", "Tissue").
#' @param sample_id unique identifier for each sample. For example "uniq_id" in the
#' RNAseqER example dataset
#' @param cols_milo_design names of columns in metadata file which are used to
#' construct the milo design matrix. Both technical and biological sources of
#' variance may be included for example c("uniq_id","Tissue","gender","AgeGroup",
#' "caseNO"). All variables of interest that will be tested for differential
#' abundance have to be added.
#' @param milo_test_fact Factors that will be used for both Milo differential
#' abundance testing and for differential gene expression analysis with
#' condition. For example c("Tissue", "gender", "AgeGroup"). Note that variables
#' here have to be an element of cols_milo_design
#' @param nfeatures number of variable features. Default is 2000.
#' @param weight_factor This is a factor used to calculate the maximum purity to
#' determine a cluster resolution that separates physically distinct clusters
#' in the umap space. The mean cluster purity is multiplied by the weight factor
#' (default is 22) and then divided by the numbers of clusters. This favors
#' clusters with a large purity over those with
#' @param pure_thres
#' @param second_thres
#'
#' @return
#' @export
#'
#' @examples
quick_RNAseqEr <- function(seur_obj,

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
                           plot_height_ann = 6,
                           dge_present_ann = FALSE,
                           dge_ann = TRUE,
                           proportion = 3,

                           #cluster qc
                           cluster_qc_vars,

                           # Milo
                           sample_id,
                           cols_milo_design,
                           milo_test_fact
                           ){
  save_dir_save <- save_dir

  #######
  ####### Set up writing to protocol
  prot_dir <- paste0(save_dir, "/outs/protocol/text_file")
  if(dir.exists(prot_dir) == FALSE){
    dir.create(prot_dir, recursive = TRUE)
  }

  #fileConn <- file(paste0(prot_dir, "/protocol.txt"))
  prot_file <- paste0(prot_dir, "/protocol.txt")

  ######
  ######



  # Perform standard Seurat processing

  seur_obj <- seurat_proc(seur_obj,
                          nfeatures = nfeatures,
                          n_pcs = n_pcs,
                          res = res,
                          select_genes = select_genes,
                          elbow_dims = elbow_dims,
                          tsne = tsne)

  ######
  ###### writing to protocol
  output_text <- paste0("Seurat processing was performed for the cell lineage ",
                        dir_lab,
                        ". The top ",
                        nfeatures,
                        " genes were chosen as variable genes with default
                        Seurat parameters for FindVariableFeatures(). ",
                        n_pcs,
                        " principal components were used for linear and non-linear
                        dimensional reduction. Check the elbow plot saved in the
                        directory outs/*cell_lineage*/plots/elbow_plot/elbow.pdf
                        whether that choice was appropriate and if not re-run with
                        a new chosen n_pcs argument. Following resolutions were
                        used for Louvain clustering: ",
                        paste(shQuote(res, type = "cmd"),
                              collapse = ", "), ".")

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- gsub("\"","",output_text_styled)
  output_text_styled <- str_squish(output_text_styled)

  #writeLines(output_text_styled, fileConn)
  sink(prot_file)
  cat(output_text_styled)
  sink()


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


  #####
  ##### write to protocol

  output_text <- paste0("The ", use_reduction, " dimensional reduction was used
                        to plot DimPlots at all found clustering resolutions that have
                        the name prefix ",
                        col_pattern,
                        ". Those plots were saved at ", save_dir, "/outs/",
                        dir_lab, "/plots/resolution_plots.")

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####

  # Calculate and plot cluster purity measures


  clu_pure(seur_obj,
           col_pattern = col_pattern,
           plot_cols = plot_cols,
           clust_lab = clust_lab_sil,
           label_size = label_size,
           save_dir = save_dir, #clu_pur_dir,
           width=7,
           height=5)



  #read in purity measures to see which cluster resolution is the one with the
  #largest number of clusters and largest purity measure

  keep_res <- max_pure(save_dir = save_dir,
                       weight_factor = weight_factor,
                       pure_thres = pure_thres,
                       second_thres = second_thres)

  Idents(seur_obj) <- keep_res

  #####
  ##### write to protocol

  output_text <- paste0("Cluster purity measures were calculated for each clustering
                        resolution using the function neighborPurity() from the bluster
                        library using their standard parameters. Output plots and
                        tables are saved at ", save_dir, "/outs/", dir_lab, "/plots/cluster_purity_plots
                        and ", save_dir, "outs/", dir_lab, "/tables/cluster_purity_data respecitively.
                        The RNAseqEr function max_pure() was used to select the most
                        appropriate clustering resolution that should distinguish
                        cell lineages (physically separated clusters in the dimensional
                        reduced space). If the result does not look like expected,
                        following parameters may be changed: weight_factor, pure_thres
                        and second_thres. The resolution being used below is: ",
                        keep_res)

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####


  ## Annotation
  seur_obj <- annotate_seqEr(seur_obj,
                             dir_lab = dir_lab,
                             ident_use = keep_res,
                             tissue_ref = tissue_ref_annotation,
                             test_use = test_use,
                             min_pct = min_pct_ann,
                             logfc_threshold = logfc_threshold_ann,
                             assay_use = assay_use,
                             save_dir = save_dir,
                             custom_ref_genes = custom_ref_genes_ann,
                             custom_gene_list,
                             customclassif = customclassif,
                             plot_reduction = plot_reduction,
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

  ####
  #### write to protocol
  output_text <- paste0("Cell type annotation was performed using the R library
                        ScType while setting the tissue reference to ", tissue_ref_annotation,
                        " . The ScType algorithm was extended to re-run it
                        also on differentially expressed genes between clusters which takes longer
                        but has improved the annotation of clusters that were previously
                        annotated as 'Unknown'. If the differential gene expression
                        analysis added valuable information, the ScType_annotation will be
                        different from the RNASeqEr_annotation. Otherwise, it they will be equal.
                        Plots showing both the ScType and RNAseqEer annotations can be found in ",
                        save_dir, "/outs/", dir_lab, "/plots/ScType_Annotated_plot.
                        The differential gene expression analysis will also be used
                        below to show celltype segregation based on gene expression. Following
                        parameters were used for the differnetial gene expression analysis: Seurat's
                        FindAllMarkers() with min.pct = ", min_pct_ann, ", test.use ", test_use,
                        ", assay = ", assay_use_ann, ", logfc.threshold = ",
                        logfc_threshold_ann,
                        ".")

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####


  #run differential gene expression for chosen clustering resolution
  Idents(seur_obj) <- "RNAseqEr_annotation"

  dge_dir <- paste0(save_dir, "/outs/", dir_lab, "/tables/DGE/broad_celltype_markers")

  clu_mark <- gen_mark_list(file_dir = dge_dir)

  int_genes <- unique(clu_mark$gene)

  df_seqEr <- heatmap_seqEr(seur_obj,
                            use_resol = keep_res,
                            dir_lab = dir_lab,
                            col_names = "RNAseqEr_annotation",
                            int_genes = int_genes,
                            save_dir = save_dir)


  ####
  #### write to protocol
  output_text <- paste0("The differential gene expression analysis used for the
                        RNAseqEr annotation at the resolution ", keep_res, " was
                        used to generate a list of top differentally expressed genes
                        between cell lineages which is used for plotting a heatmap
                        which is saved in outs/", dir_lab, "/plots/heatmaps.")

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####

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

  ####
  #### write to protocol
  output_text <- paste0("Cluster quality control was performed by using the
                        RNAseqEer function cluster_qc() which looks for
                        the proportion of cells/ nuclei per cluster derived from
                        each one of the following variables of interest: ",
                        paste(shQuote(cluster_qc_vars, type = "cmd"),
                              collapse = ", "), " using the cluster labels in ",
                        cluster_col,
                        ". This is done to look for technical
                        not biological variation as clusters that fail quality
                        control thresholds are being removed from the dataset
                        using the function remove_clu() from the RNAseqEr library.
                        The default settings for cluster_qc() are currently used
                        which fail a cluster if not at least 70% of each variable
                        provided in cluster_qc_vars contributes to each cluster. If
                        this filtering seems to be not right for the data (check
                        bargraphs at outs/*cell_type*/plots/cluster_qc, DimPlots
                        at outs/*cell_type*/plots/cluster_qc/ DimPlots and data at
                        outs/*cell_type*/tables/cluster_qc), get in touch
                        (lseeker@stanford.edu). The data is used to filter the data
                        and only those that have 'passed' the QC are retained."
                        )

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- gsub("\"","",output_text_styled)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####


  # Test if abundancy differs with variables of technical or biological interest
  milo_obj <- abundance_test(seur_obj = seur_obj,
                             cluster_col = "RNAseqEr_annotation",
                             sample_id,
                             cols_interest = cols_milo_design,
                             test_factors = milo_test_fact,
                             dir_lab = dir_lab)


  ####
  #### write to protocol
  output_text <- paste0("Subsequently Milo (Dann et al. 2022) is used to
                        test for variation in abundance of cells/ nuclei using the
                        RNAseqEr function abundance_test() with default settings
                        and following variables of interest: ",
                        paste(shQuote(milo_test_fact, type = "cmd"),
                              collapse = ", "),
                        ". All variables of interest which may be biological or
                        technical should be added to 'cols_milo_design'. In this
                        case those were: ",
                        paste(shQuote(cols_milo_design, type = "cmd"),
                              collapse = ", "),
                        ". Following sample ID was used: ",
                        sample_id,
                        " to identify single samples individual samples."

  )

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- gsub("\"","",output_text_styled)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####




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

  ####
  #### write to protocol
  output_text <- paste0("The RNAseqEr function find_cond_markers() is used to
                        look for variation in gene expression with following biological
                        variables of interest: ",
                        paste(shQuote(milo_test_fact, type = "cmd"),
                              collapse = ", "),
                        ". The FindAllMarkers() function from the Seurat libary is
                        used with paremeters default paremeters and following
                        parameters that can be adapted and are currently set to
                        the values below:  logfc.threshold = ",
                        logfc_threshold,
                        ", only.pos = ",
                        only_pos,
                        ", min.pct = ",
                        min_pct,
                        ", test.use = ",
                        test_use,
                        ". The differential gene expression analysis is perfommed
                        once across the entire dataset
                        (outs/*cell_linage*/tables/condition_mark/RNAseqEr_annotation/overall)
                        and once clusterwise
                        (outs/*cell_linage*/tables/condition_mark/RNAseqEr_annotation/clusterwise),
                        using cluster labels in the column ",
                        cluster_id,
                        " which at this point should correspond to cell lineage labels.
                        The differential gene expression result of the overall test
                        (not clusterwise) is filtered based on the proportion of cells within the
                        cluster of interest expressing the gene being pct1 > ",
                        fil_pct_1,
                        " and the proportion of cells outside the cluster of interest
                        expressing the gene being pct2 < ",
                        fil_pct_2,
                        ". The filtered data is saved separately at using a 'fil_'
                        instead of an 'all_' prefix for the file name.
                        The differential gene expression results are used by the
                        RNAseqEr function gen_mark_list() for each biological variable
                        of interest to filter for the top 10 genes per cluster. Those
                        genes are plotted as violin and feature plots using Seurat.
                        Those plots can be found here:
                        outs/*cell_type*/plots/Condition/*variable_name*_mark/violin_plots/RNAseqER_annotation and
                        outs/*cell_type*/plots/Condition/*variable_name*_mark/FeaturePlots."
                        )

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- gsub("\"","",output_text_styled)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####

  # Now that we have a resolution for the dataset that enabled a preliminary annotatoion,
  # let's subset for the main clusters (as in cell lineages) and look for
  # finer clusters that make sense biologically


  subset_RNAseqEr(seur_obj,
                  subset_column = "RNAseqEr_annotation",
                  save_dir = save_dir)






  # process cell lineage datasets

  save_path <- paste0(save_dir, "/outs/data")
  cl_dat_name <- list.files(save_path)

  for(y in 1: length(cl_dat_name)){
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



    df_seqEr <- heatmap_seqEr(curr_srt,
                              dir_lab = cell_lineage_name,
                              use_resol = use_resol,
                              col_names = keep_res_ul,
                              int_genes = int_genes,
                              save_dir = save_dir,
                              max_diff_threshold = max_diff_threshold,
                              mean_diff_thres = mean_diff_thres)

    df_seqEr$res <- as.numeric(sapply(strsplit(df_seqEr$resolution, col_pattern), "[", 2))

    sub_df <- subset(df_seqEr, df_seqEr$cluster_sim_score == min(df_seqEr$cluster_sim_score))
    sub_df <- subset(sub_df, sub_df$res == max(sub_df$res))

    keep_res <- sub_df$resolution


    Idents(curr_srt) <- keep_res

    curr_srt$cluster_id <- paste0(cell_lineage_name, "_", curr_srt@meta.data[[keep_res]])

    annot_dim <- DimPlot(curr_srt, group.by = "cluster_id", label = TRUE, cols = colour_palette())

    dimplot_dir <- paste0(save_dir,
                          "/outs/",
                          cell_lineage_name,
                          "/plots/dim_plot_annotated")
    if(dir.exists(dimplot_dir) == FALSE){
      dir.create(dimplot_dir, recursive = TRUE)
    }

    pdf(paste0(dimplot_dir, "/anoot_dim.pdf"),
        height=plotheight, width = plotwidth)
    print(annot_dim)
    dev.off()

    print(annot_dim)


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
        #clu_levels <- levels(as.factor(mark$cluster))

        #for(u in 1: length(clu_levels)){
        #  subs_dat <- subset(mark, mark$cluster == clu_levels[u])
         # curr_level <- clu_levels[u]

         # go_results <- perform_go(curr_srt,
         #                          gene_list = mark,
         #                          min_log2FC = 0.25,
         #                          reverse = FALSE,
         #                          translate_gene_id_from = "SYMBOL",
         #                          translate_gene_id_to = "ENTREZID",
          #                         org_use = "org.Hs.eg.db",
          #                         ontology = "BP",
          #                         pvalue_cutoff = 0.05,
          #                         qvalue_cutoff = 0.05,
          #                         read_able = TRUE)

          #go_plot <- dotplot(go_results)

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

  ####
  #### write to protocol
  output_text <- paste0("The dataset was divided into smaller cell lineage datasets
                        and saved as separate .RDS objects in outs/data.
                        Each of those cell lineage datasets was then individually
                        processed. First, the RNAseqEr function seurat_proc()
                        was used with the same parameters as for the complete
                        dataset. The selection of the suitable cluster resolution is
                        performed differentty for cell lineage datasets than for the
                        complete dataset. Instead of relying on purity measues,
                        A cluster resolution within cell lineages is searched that
                        allows for the cluster segregation based on differentially
                        expressed genes. The RNAseqEr function select_res() is
                        used to select resolutions that are being used for differential
                        gene expression including a pairwise comparison of all
                        cluster levels. This step takes significant time therefore the
                        of overclustering resolutions with an excessive number of clusters
                        should be avoided. All the resulting data from the differntial
                        gene expression analyses (saved at
                        outs/*cell_linage*/tables/cluster_marker) is processed and
                        filtered for the top genes separating clusters (default is
                        top 10 genes per cluster). The results are plotted in
                        heatmaps that are saved at
                        outs/*cell_linage*/heatmaps. Those heatmaps can be visually
                        inpected for cluster segregation which is also done using the
                        function heatmap_seqEr(). The cluster quality control is performed in paralell to
                        the cluster quality control for the entire dataset with the
                        dfference that this time flagged clusters are not automatically
                        removes as at this level they are more likely to be affected by biological
                        variance. Subsequently, Milo (Dann et al. 2022) is used for
                        a differential abundance test, using the same factors for
                        building the design matrix and for testing as for the
                        complete dataset: ",
                        paste(shQuote(milo_test_fact, type = "cmd"),
                              collapse = ", "),
                        " Also in paralell to the complete dataset, a differential gene
                        expression analysis with those conditions is performed once across the
                        entire cell-linage dataset, once cluster-wise.)"
  )

  output_text_styled <- gsub("\n","",output_text)
  output_text_styled <- gsub("\"","",output_text_styled)
  output_text_styled <- str_squish(output_text_styled)

  sink(prot_file, append = TRUE)
  cat(output_text_styled)
  sink()

  ####
  ####



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
                                    comp_mode = comp_modes[v],
                                    conditions = milo_test_fact,
                                    cluster_label = "cell_lineage")
    volcano_seqEr(curr_mark,clu_col = "cluster_id",
                              condition_label = "condition",
                              mode = comp_modes[v],
                              save_dir = getwd(),
                              conditions = milo_test_fact,
                              plotheight = 6,
                              plotwidth = 8,
                              save_plots = TRUE)


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

  close(fileConn)
  return(seur_obj)

 }






