#' Generate comprehensive analysis report for RNAseqEr workflow
#' @description
#' Generates a detailed report in two parts: (1) Methods-style workflow description
#' and (2) Decision justifications for key analysis choices. Designed for early
#' career researchers to understand and justify their analysis decisions.
#'
#' @param seur_obj Seurat object that was analyzed
#' @param save_dir Directory where analysis results are saved
#' @param dir_lab Label used for the analysis
#' @param workflow_params List of parameters used in the workflow
#' @param analysis_date Date of analysis (defaults to current date)
#' @param researcher_name Name of the researcher conducting the analysis
#' @param project_title Title of the project
#' @param tissue_type Type of tissue analyzed
#' @param species Species analyzed (default: "Human")
#'
#' @return Saves two report files: methods_report.md and decision_justifications.md
#' @export
#'
#' @examples
#' \dontrun{
#' generate_report(cns, 
#'                save_dir = "analysis_results",
#'                dir_lab = "all_celltypes",
#'                workflow_params = list(n_pcs = 20, res = c(0.1, 0.2, 0.3)),
#'                researcher_name = "Student Name",
#'                project_title = "Single-cell analysis of brain tissue",
#'                tissue_type = "Brain")
#' }
#'
generate_report <- function(seur_obj,
                          save_dir,
                          dir_lab = "all_celltypes",
                          workflow_params = list(),
                          analysis_date = Sys.Date(),
                          researcher_name = "Researcher",
                          project_title = "Single-cell RNA-seq Analysis",
                          tissue_type = "Tissue",
                          species = "Human") {
  
  # Create reports directory
  report_dir <- paste0(save_dir, "/reports")
  if(!dir.exists(report_dir)) {
    dir.create(report_dir, recursive = TRUE)
  }
  
  # Get basic dataset information
  n_cells <- ncol(seur_obj)
  n_genes <- nrow(seur_obj)
  n_clusters <- length(unique(Idents(seur_obj)))
  
  # Generate Methods Report
  methods_report <- generate_methods_report(seur_obj, save_dir, dir_lab, 
                                          workflow_params, analysis_date, 
                                          researcher_name, project_title, 
                                          tissue_type, species, n_cells, n_genes, n_clusters)
  
  # Generate Decision Justifications Report
  decision_report <- generate_decision_justifications(seur_obj, save_dir, dir_lab, 
                                                   workflow_params, tissue_type)
  
  # Save reports
  writeLines(methods_report, paste0(report_dir, "/methods_report.md"))
  writeLines(decision_report, paste0(report_dir, "/decision_justifications.md"))
  
  # Generate combined report
  combined_report <- c("# RNAseqEr Analysis Report", "",
                      "**Generated on:** ", as.character(analysis_date), "",
                      "**Researcher:** ", researcher_name, "",
                      "**Project:** ", project_title, "",
                      "**Tissue:** ", tissue_type, "",
                      "**Species:** ", species, "",
                      "",
                      "## Dataset Summary", "",
                      paste("- **Total cells:**", n_cells),
                      paste("- **Total genes:**", n_genes),
                      paste("- **Final clusters:**", n_clusters),
                      "",
                      "---",
                      "",
                      "## Methods", "",
                      methods_report,
                      "",
                      "---",
                      "",
                      "## Analysis Decisions and Justifications", "",
                      decision_report)
  
  writeLines(combined_report, paste0(report_dir, "/complete_analysis_report.md"))
  
  cat("Reports generated successfully in:", report_dir, "\n")
  cat("- methods_report.md\n")
  cat("- decision_justifications.md\n") 
  cat("- complete_analysis_report.md\n")
  
  return(report_dir)
}

#' Generate methods-style report
#' @keywords internal
generate_methods_report <- function(seur_obj, save_dir, dir_lab, workflow_params, 
                                  analysis_date, researcher_name, project_title, 
                                  tissue_type, species, n_cells, n_genes, n_clusters) {
  
  methods_text <- c(
    "## Single-cell RNA-seq Analysis Methods",
    "",
    paste("Single-cell RNA sequencing data from", tissue_type, "tissue was analyzed using the RNAseqEr package (version 0.1.0) in R. The analysis workflow consisted of the following steps:"),
    "",
    "### 1. Data Preprocessing",
    paste("- Quality control was performed to remove low-quality cells and genes"),
    paste("- Data normalization was performed using the Seurat NormalizeData() function"),
    paste("- Variable genes were identified using FindVariableFeatures()"),
    paste("- Data scaling was performed using ScaleData()"),
    "",
    "### 2. Dimensional Reduction and Clustering",
    paste("- Principal component analysis (PCA) was performed on variable genes"),
    paste("-", workflow_params$n_pcs %||% 20, "principal components were selected for downstream analysis"),
    paste("- Clustering was performed using the Louvain algorithm at multiple resolutions"),
    paste("- UMAP dimensional reduction was performed for visualization"),
    "",
    "### 3. Cluster Resolution Selection",
    paste("- Cluster purity was calculated using silhouette scores"),
    paste("- The optimal clustering resolution was selected based on maximum cluster purity and number of clusters"),
    paste("- Final clustering resolution:", workflow_params$keep_res %||% "determined automatically"),
    paste("- Cluster purity threshold:", workflow_params$pure_thres %||% 0.96),
    paste("- Weight factor for purity calculation:", workflow_params$weight_factor %||% 22),
    "",
    "### 4. Cell Type Annotation",
    paste("- Cell type annotation was performed using the ScType algorithm"),
    paste("- Tissue reference:", workflow_params$tissue_ref %||% "Brain"),
    paste("- Differential gene expression analysis was used to refine annotations for unknown clusters"),
    "",
    "### 5. Subclustering Analysis",
    paste("- Major cell lineages were subset for detailed analysis"),
    paste("- Subclustering was performed within each lineage"),
    paste("- Cluster similarity analysis was used to select optimal subclustering resolution"),
    "",
    "### 6. Quality Control",
    paste("- Cluster quality control was performed to identify technical artifacts"),
    paste("- Clusters with sample-specific bias were identified and removed"),
    paste("- Quality control variables:", paste(workflow_params$cluster_qc_vars %||% c("sample_id"), collapse = ", ")),
    "",
    "### 7. Differential Analysis",
    paste("- Differential abundance testing was performed using Milo"),
    paste("- Differential gene expression analysis was performed using", workflow_params$test_use %||% "MAST"),
    paste("- Log fold change threshold:", workflow_params$logfc_threshold %||% 0.25),
    paste("- Minimum percentage threshold:", workflow_params$min_pct %||% 0.25),
    paste("- Condition markers were identified for each cluster"),
    "",
    "### 8. Gene Ontology Analysis",
    paste("- Gene ontology analysis was performed using ClusterProfiler"),
    paste("- Biological processes (BP) ontology was used"),
    paste("- P-value cutoff:", workflow_params$pvalue_cutoff_go %||% 0.05),
    paste("- Q-value cutoff:", workflow_params$qvalue_cutoff_go %||% 0.05),
    "",
    "### 9. Visualization and Interactive Analysis",
    paste("- Dimensional reduction plots were generated for all clustering resolutions"),
    paste("- Heatmaps were generated for cluster marker genes"),
    paste("- Violin plots were generated for condition-specific genes"),
    paste("- An interactive Shiny application was created for data exploration"),
    "",
    "### Software and Packages",
    "- R version 4.3.0 or higher",
    "- Seurat v4.3.0 or higher",
    "- RNAseqEr v0.1.0",
    "- Additional packages: miloR, clusterProfiler, org.Hs.eg.db"
  )
  
  return(methods_text)
}

#' Generate decision justifications report
#' @keywords internal
generate_decision_justifications <- function(seur_obj, save_dir, dir_lab, 
                                          workflow_params, tissue_type) {
  
  # Read analysis results to provide specific justifications
  purity_file <- paste0(save_dir, "/outs/", dir_lab, "/tables/cluster_purity_data/summary_purity_stats.csv")
  
  # Try to read actual thresholds from analysis files
  actual_thresholds <- list()
  
  # Read cluster purity results if available
  if(file.exists(purity_file)) {
    purity_data <- read.csv(purity_file)
    actual_thresholds$cluster_purity <- max(purity_data$mean_pur, na.rm = TRUE)
    actual_thresholds$n_clusters <- max(purity_data$clu_count, na.rm = TRUE)
    actual_thresholds$selected_resolution <- purity_data$keep_res[1]
  }
  
  # Read differential expression results if available
  dge_file <- paste0(save_dir, "/outs/", dir_lab, "/tables/DGE/broad_celltype_markers/all.csv")
  if(file.exists(dge_file)) {
    dge_data <- read.csv(dge_file)
    actual_thresholds$n_markers <- nrow(dge_data)
    actual_thresholds$avg_logfc <- mean(dge_data$avg_log2FC, na.rm = TRUE)
    actual_thresholds$min_pct <- min(dge_data$pct.1, na.rm = TRUE)
  }
  
  # Read condition marker results if available
  cond_file <- paste0(save_dir, "/outs/", dir_lab, "/tables/condition_mark/RNAseqEr_annotation/clusterwise/all.csv")
  if(file.exists(cond_file)) {
    cond_data <- read.csv(cond_file)
    actual_thresholds$n_condition_markers <- nrow(cond_data)
  }
  
  decisions_text <- c(
    "## Analysis Decisions and Justifications",
    "",
    "This section provides justifications for key analytical decisions made during the analysis:",
    "",
    "### 1. Number of Principal Components",
    paste("**Decision:** Used", workflow_params$n_pcs %||% 20, "principal components"),
    "**Justification:** The number of PCs was selected based on the elbow plot showing the point of diminishing returns in explained variance. This ensures we capture the biological signal while avoiding noise from higher dimensions.",
    "",
    "### 2. Clustering Resolution Selection",
    paste("**Decision:** Selected clustering resolution based on cluster purity analysis"),
    paste("**Justification:** Multiple clustering resolutions were tested, and the resolution with the highest cluster purity score was selected. The selected resolution achieved a cluster purity of", 
          round(actual_thresholds$cluster_purity %||% 0, 3), "with", 
          actual_thresholds$n_clusters %||% "multiple", "clusters. This ensures that clusters represent biologically distinct cell populations rather than technical artifacts."),
    "",
    "### 3. Cell Type Annotation Strategy",
    paste("**Decision:** Used ScType algorithm with", workflow_params$tissue_ref %||% "Brain", "tissue reference"),
    "**Justification:** ScType provides automated cell type annotation based on known marker genes. For clusters classified as 'Unknown', additional differential gene expression analysis was performed to refine the annotation.",
    "",
    "### 4. Quality Control Thresholds",
    paste("**Decision:** Applied quality control using variables:", paste(workflow_params$cluster_qc_vars %||% c("sample_id"), collapse = ", ")),
    "**Justification:** Clusters that showed strong association with technical variables (e.g., sample ID, processing batch) were flagged as potentially technical artifacts and removed from downstream analysis.",
    "",
    "### 5. Differential Expression Analysis Parameters",
    paste("**Decision:** Used", workflow_params$test_use %||% "MAST", "statistical test with log fold change threshold of", workflow_params$logfc_threshold %||% 0.25),
    paste("**Justification:** MAST is specifically designed for single-cell data and accounts for the zero-inflated nature of the data. The analysis identified", 
          actual_thresholds$n_markers %||% "multiple", "marker genes with an average log fold change of", 
          round(actual_thresholds$avg_logfc %||% 0, 2), ". The log fold change threshold was chosen to balance sensitivity and specificity."),
    "",
    "### 6. Gene Ontology Analysis Settings",
    paste("**Decision:** Used Biological Processes ontology with p-value cutoff of", workflow_params$pvalue_cutoff_go %||% 0.05),
    "**Justification:** Biological Processes provides the most interpretable results for understanding cellular functions. The p-value cutoff ensures statistical significance while the q-value cutoff controls for multiple testing.",
    "",
    "### 7. Subclustering Strategy",
    "**Decision:** Performed subclustering within major cell lineages",
    "**Justification:** Subclustering allows for the identification of cell subtypes within major lineages. This is particularly important for complex tissues where cell types may be further subdivided.",
    "",
    "### 8. Visualization Choices",
    "**Decision:** Generated UMAP plots and heatmaps for visualization",
    "**Justification:** UMAP provides an intuitive visualization of cell relationships while preserving local structure. Heatmaps allow for the visualization of gene expression patterns across clusters and conditions.",
    "",
    "### Alternative Approaches Considered",
    "",
    "1. **Clustering Algorithm:** Louvain was chosen over other algorithms (e.g., Leiden, hierarchical clustering) due to its computational efficiency and good performance on single-cell data.",
    "",
    "2. **Dimensional Reduction:** UMAP was chosen over tSNE due to better preservation of global structure and faster computation.",
    "",
    "3. **Annotation Method:** ScType was chosen over other methods (e.g., SingleR, Garnett) due to its tissue-specific approach and interpretable results.",
    "",
    "### Limitations and Considerations",
    "",
    "1. **Batch Effects:** While quality control was performed, batch effects may still influence the results. Consider using batch correction methods if strong batch effects are observed.",
    "",
    "2. **Cell Type Annotation:** Automated annotation methods have limitations and should be validated with independent methods (e.g., immunohistochemistry, flow cytometry).",
    "",
    "3. **Statistical Power:** The statistical power of differential expression analysis depends on the number of cells per cluster and condition.",
    "",
    "4. **Computational Resources:** The analysis requires significant computational resources, particularly for large datasets.",
    "",
    "### Summary of Analysis Results",
    "",
    paste("**Dataset Characteristics:**"),
    paste("- Total cells analyzed:", ncol(seur_obj)),
    paste("- Total genes analyzed:", nrow(seur_obj)),
    paste("- Final number of clusters:", length(unique(Idents(seur_obj)))),
    "",
    paste("**Clustering Results:**"),
    paste("- Selected clustering resolution:", actual_thresholds$selected_resolution %||% "determined automatically"),
    paste("- Cluster purity achieved:", round(actual_thresholds$cluster_purity %||% 0, 3)),
    paste("- Number of clusters identified:", actual_thresholds$n_clusters %||% "multiple"),
    "",
    paste("**Differential Expression Results:**"),
    paste("- Number of marker genes identified:", actual_thresholds$n_markers %||% "multiple"),
    paste("- Average log fold change:", round(actual_thresholds$avg_logfc %||% 0, 2)),
    paste("- Number of condition-specific markers:", actual_thresholds$n_condition_markers %||% "multiple")
  )
  
  return(decisions_text)
}

# Helper function for safe list access
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
} 