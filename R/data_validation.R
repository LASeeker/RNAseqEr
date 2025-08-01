#' Check cluster sizes and validate for differential expression analysis
#' @description
#' Validates that clusters have sufficient cells for reliable differential expression analysis
#'
#' @param seur_obj Seurat object
#' @param cluster_col Column name containing cluster assignments
#' @param min_cells_per_cluster Minimum number of cells required per cluster (default: 3)
#' @param min_cells_for_comparison Minimum cells needed for pairwise comparisons (default: 5)
#'
#' @return List with validation results and recommendations
#' @export
#'
#' @examples
#' validation_result <- validate_cluster_sizes(seur_obj, "RNA_snn_res.0.8")
#'
validate_cluster_sizes <- function(seur_obj, 
                                  cluster_col, 
                                  min_cells_per_cluster = 3,
                                  min_cells_for_comparison = 5) {
  
  # Get cluster sizes
  cluster_sizes <- table(seur_obj@meta.data[[cluster_col]])
  
  # Identify problematic clusters
  small_clusters <- cluster_sizes < min_cells_per_cluster
  very_small_clusters <- cluster_sizes < min_cells_for_comparison
  
  # Calculate total cells
  total_cells <- sum(cluster_sizes)
  
  # Check if we have enough clusters for comparison
  n_clusters <- length(cluster_sizes)
  n_valid_clusters <- sum(cluster_sizes >= min_cells_per_cluster)
  
  # Determine if analysis is feasible
  feasible <- n_valid_clusters >= 2
  
  # Create recommendations
  recommendations <- list()
  
  if(!feasible) {
    recommendations$action <- "skip_analysis"
    recommendations$reason <- "Insufficient clusters with minimum cells"
  } else if(sum(small_clusters) > 0) {
    recommendations$action <- "filter_clusters"
    recommendations$small_clusters <- names(cluster_sizes[small_clusters])
    recommendations$reason <- "Some clusters have too few cells"
  } else {
    recommendations$action <- "proceed"
    recommendations$reason <- "All clusters have sufficient cells"
  }
  
  return(list(
    cluster_sizes = cluster_sizes,
    total_cells = total_cells,
    n_clusters = n_clusters,
    n_valid_clusters = n_valid_clusters,
    small_clusters = names(cluster_sizes[small_clusters]),
    very_small_clusters = names(cluster_sizes[very_small_clusters]),
    feasible = feasible,
    recommendations = recommendations
  ))
}

#' Filter Seurat object to remove small clusters
#' @description
#' Removes clusters with insufficient cells for analysis
#'
#' @param seur_obj Seurat object
#' @param cluster_col Column name containing cluster assignments
#' @param min_cells_per_cluster Minimum number of cells required per cluster
#' @param verbose Whether to print information about filtering
#'
#' @return Filtered Seurat object
#' @export
#'
#' @examples
#' filtered_obj <- filter_small_clusters(seur_obj, "RNA_snn_res.0.8", min_cells_per_cluster = 3)
#'
filter_small_clusters <- function(seur_obj, 
                                 cluster_col, 
                                 min_cells_per_cluster = 3,
                                 verbose = TRUE) {
  
  # Validate cluster sizes
  validation <- validate_cluster_sizes(seur_obj, cluster_col, min_cells_per_cluster)
  
  if(validation$recommendations$action == "proceed") {
    if(verbose) {
      cat("✅ All clusters have sufficient cells. No filtering needed.\n")
    }
    return(seur_obj)
  }
  
  if(validation$recommendations$action == "skip_analysis") {
    if(verbose) {
      cat("❌ Insufficient clusters for analysis. Consider:\n")
      cat("   - Lowering clustering resolution\n")
      cat("   - Increasing minimum cells per cluster\n")
      cat("   - Using different clustering parameters\n")
    }
    return(NULL)
  }
  
  # Filter out small clusters
  small_clusters <- validation$recommendations$small_clusters
  cells_to_keep <- !(seur_obj@meta.data[[cluster_col]] %in% small_clusters)
  
  if(verbose) {
    cat("⚠️  Filtering out small clusters:", paste(small_clusters, collapse = ", "), "\n")
    cat("   Cells removed:", sum(!cells_to_keep), "\n")
    cat("   Cells remaining:", sum(cells_to_keep), "\n")
  }
  
  # Subset the object
  filtered_obj <- seur_obj[, cells_to_keep]
  
  # Recalculate cluster assignments if needed
  if(length(unique(filtered_obj@meta.data[[cluster_col]])) < 2) {
    if(verbose) {
      cat("⚠️  After filtering, only one cluster remains. Analysis may not be meaningful.\n")
    }
  }
  
  return(filtered_obj)
}

#' Safe differential expression analysis with validation
#' @description
#' Performs differential expression analysis with safeguards for small datasets
#'
#' @param seur_obj Seurat object
#' @param cluster_col Column name containing cluster assignments
#' @param min_cells_per_cluster Minimum cells per cluster (default: 3)
#' @param min_cells_for_comparison Minimum cells for pairwise comparisons (default: 5)
#' @param verbose Whether to print progress information
#'
#' @return List with analysis results and warnings
#' @export
#'
#' @examples
#' safe_dge_result <- safe_differential_expression(seur_obj, "RNA_snn_res.0.8")
#'
safe_differential_expression <- function(seur_obj, 
                                       cluster_col,
                                       min_cells_per_cluster = 3,
                                       min_cells_for_comparison = 5,
                                       verbose = TRUE) {
  
  # Validate cluster sizes
  validation <- validate_cluster_sizes(seur_obj, cluster_col, 
                                      min_cells_per_cluster, min_cells_for_comparison)
  
  if(verbose) {
    cat("=== Data Validation ===\n")
    cat("Total cells:", validation$total_cells, "\n")
    cat("Number of clusters:", validation$n_clusters, "\n")
    cat("Valid clusters:", validation$n_valid_clusters, "\n")
    
    if(length(validation$small_clusters) > 0) {
      cat("Small clusters:", paste(validation$small_clusters, collapse = ", "), "\n")
    }
  }
  
  # Check if analysis is feasible
  if(!validation$feasible) {
    warning("Differential expression analysis skipped: insufficient cells/clusters")
    return(list(
      success = FALSE,
      reason = "Insufficient cells or clusters for analysis",
      validation = validation
    ))
  }
  
  # Filter if needed
  if(validation$recommendations$action == "filter_clusters") {
    if(verbose) {
      cat("Filtering small clusters...\n")
    }
    seur_obj <- filter_small_clusters(seur_obj, cluster_col, min_cells_per_cluster, verbose)
    
    if(is.null(seur_obj)) {
      return(list(
        success = FALSE,
        reason = "No valid clusters after filtering",
        validation = validation
      ))
    }
  }
  
  # Perform analysis with error handling
  tryCatch({
    if(verbose) {
      cat("Performing differential expression analysis...\n")
    }
    
    # Set identities
    Idents(seur_obj) <- cluster_col
    
    # Perform FindAllMarkers with error handling
    all_markers <- FindAllMarkers(seur_obj,
                                  only.pos = TRUE,
                                  min.pct = 0.25,
                                  logfc.threshold = 0.25,
                                  test.use = "MAST",
                                  verbose = FALSE)
    
    if(verbose) {
      cat("✅ Differential expression analysis completed successfully\n")
      cat("   Markers found:", nrow(all_markers), "\n")
    }
    
    return(list(
      success = TRUE,
      markers = all_markers,
      validation = validation,
      filtered_object = seur_obj
    ))
    
  }, error = function(e) {
    warning("Differential expression analysis failed: ", e$message)
    return(list(
      success = FALSE,
      reason = e$message,
      validation = validation
    ))
  })
}

#' Check dataset size and provide recommendations
#' @description
#' Analyzes dataset size and provides recommendations for analysis parameters
#'
#' @param seur_obj Seurat object
#' @param verbose Whether to print detailed recommendations
#'
#' @return List with size analysis and recommendations
#' @export
#'
#' @examples
#' size_analysis <- analyze_dataset_size(seur_obj)
#'
analyze_dataset_size <- function(seur_obj, verbose = TRUE) {
  
  total_cells <- ncol(seur_obj)
  
  # Size categories
  if(total_cells < 100) {
    size_category <- "very_small"
    recommendations <- list(
      clustering_resolution = "high (0.1-0.5)",
      min_cells_per_cluster = 3,
      min_cells_for_comparison = 5,
      statistical_test = "wilcox",
      warning = "Very small dataset - results may not be reliable"
    )
  } else if(total_cells < 1000) {
    size_category <- "small"
    recommendations <- list(
      clustering_resolution = "medium (0.3-0.8)",
      min_cells_per_cluster = 5,
      min_cells_for_comparison = 10,
      statistical_test = "MAST",
      warning = "Small dataset - consider conservative thresholds"
    )
  } else if(total_cells < 10000) {
    size_category <- "medium"
    recommendations <- list(
      clustering_resolution = "standard (0.1-1.5)",
      min_cells_per_cluster = 10,
      min_cells_for_comparison = 20,
      statistical_test = "MAST",
      warning = NULL
    )
  } else {
    size_category <- "large"
    recommendations <- list(
      clustering_resolution = "standard (0.1-1.5)",
      min_cells_per_cluster = 20,
      min_cells_for_comparison = 50,
      statistical_test = "MAST",
      warning = NULL
    )
  }
  
  if(verbose) {
    cat("=== Dataset Size Analysis ===\n")
    cat("Total cells:", total_cells, "\n")
    cat("Size category:", size_category, "\n")
    cat("\nRecommendations:\n")
    cat("- Clustering resolution:", recommendations$clustering_resolution, "\n")
    cat("- Min cells per cluster:", recommendations$min_cells_per_cluster, "\n")
    cat("- Min cells for comparison:", recommendations$min_cells_for_comparison, "\n")
    cat("- Statistical test:", recommendations$statistical_test, "\n")
    
    if(!is.null(recommendations$warning)) {
      cat("\n⚠️  Warning:", recommendations$warning, "\n")
    }
  }
  
  return(list(
    total_cells = total_cells,
    size_category = size_category,
    recommendations = recommendations
  ))
} 