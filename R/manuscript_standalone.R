#' Generate manuscript draft from RNAseqEr analysis results
#' @description
#' Standalone function to generate manuscript drafts from existing RNAseqEr analysis results.
#' This function can be used independently of the full quick_RNAseqEr workflow.
#'
#' @param save_dir Directory where RNAseqEr analysis results are saved
#' @param dir_lab Label used for the analysis (default: "all_celltypes")
#' @param conditions_of_interest Vector of condition variables to focus on (e.g., c("AgeGroup", "Tissue"))
#' @param tissue_type Type of tissue analyzed (default: "Brain")
#' @param species Species analyzed (default: "Human")
#' @param project_title Title for the manuscript (default: "Single-cell RNA-seq Analysis")
#' @param researcher_name Name of the researcher (default: "Researcher")
#' @param journal_target Target journal for manuscript style (default: "bioRxiv")
#' @param llm_provider LLM provider to use ("openai", "anthropic", or "local")
#' @param api_key API key for LLM provider (if NULL, will try to load from environment)
#' @param model_name Model name to use (default: "gpt-4")
#' @param max_tokens Maximum tokens for LLM response (default: 4000)
#' @param temperature Temperature for LLM generation (default: 0.7)
#' @param include_figures Whether to include figure selection (default: TRUE)
#' @param focus_on_genes Focus on specific gene categories ("condition", "cluster", or "all")
#' @param min_logfc Minimum log fold change for gene selection (default: 0.5)
#' @param max_pvalue Maximum p-value for gene selection (default: 0.05)
#' @param auto_load_env Whether to automatically load environment variables (default: TRUE)
#'
#' @return Path to the manuscript directory
#' @export
#'
#' @examples
#' \dontrun{
#' # Generate manuscript from existing analysis
#' manuscript_dir <- generate_manuscript_standalone(
#'   save_dir = "analysis_results",
#'   dir_lab = "all_celltypes",
#'   conditions_of_interest = c("AgeGroup", "Tissue"),
#'   project_title = "Age-related changes in brain tissue",
#'   researcher_name = "Dr. Smith",
#'   journal_target = "Nature Communications",
#'   llm_provider = "openai",
#'   api_key = "your-api-key-here"
#' )
#' 
#' # Use with environment variables
#' manuscript_dir <- generate_manuscript_standalone(
#'   save_dir = "analysis_results",
#'   conditions_of_interest = c("AgeGroup"),
#'   auto_load_env = TRUE  # Will load from .env file
#' )
#' }
#'
generate_manuscript_standalone <- function(save_dir,
                                         dir_lab = "all_celltypes",
                                         conditions_of_interest = NULL,
                                         tissue_type = "Brain",
                                         species = "Human",
                                         project_title = "Single-cell RNA-seq Analysis",
                                         researcher_name = "Researcher",
                                         journal_target = "bioRxiv",
                                         llm_provider = "openai",
                                         api_key = NULL,
                                         model_name = "gpt-4",
                                         max_tokens = 4000,
                                         temperature = 0.7,
                                         include_figures = TRUE,
                                         focus_on_genes = "condition",
                                         min_logfc = 0.5,
                                         max_pvalue = 0.05,
                                         auto_load_env = TRUE) {
  
  # Check if analysis results exist
  if(!dir.exists(save_dir)) {
    stop("Analysis directory not found: ", save_dir)
  }
  
  # Check for required subdirectories
  required_dirs <- c(
    paste0(save_dir, "/outs/", dir_lab, "/tables/condition_mark"),
    paste0(save_dir, "/outs/", dir_lab, "/tables/cluster_marker"),
    paste0(save_dir, "/outs/", dir_lab, "/plots")
  )
  
  missing_dirs <- required_dirs[!dir.exists(required_dirs)]
  if(length(missing_dirs) > 0) {
    warning("Some required directories are missing. Manuscript generation may be incomplete:\n",
            paste(missing_dirs, collapse = "\n"))
  }
  
  # Auto-load environment variables if requested
  if(auto_load_env && is.null(api_key)) {
    tryCatch({
      load_env_vars()
      api_key <- get_api_key(provider = llm_provider)
      if(is.null(api_key)) {
        warning("No API key found in environment. Please set api_key parameter or configure environment variables.")
      }
    }, error = function(e) {
      warning("Could not load environment variables: ", e$message)
    })
  }
  
  # Generate manuscript
  manuscript_dir <- generate_manuscript(
    save_dir = save_dir,
    dir_lab = dir_lab,
    conditions_of_interest = conditions_of_interest,
    tissue_type = tissue_type,
    species = species,
    project_title = project_title,
    researcher_name = researcher_name,
    journal_target = journal_target,
    llm_provider = llm_provider,
    api_key = api_key,
    model_name = model_name,
    max_tokens = max_tokens,
    temperature = temperature,
    include_figures = include_figures,
    focus_on_genes = focus_on_genes,
    min_logfc = min_logfc,
    max_pvalue = max_pvalue,
    auto_load_env = FALSE  # Already handled above
  )
  
  cat("✅ Manuscript draft generated successfully!\n")
  cat("📁 Manuscript files saved to:", manuscript_dir, "\n")
  cat("📄 Files created:\n")
  cat("   - manuscript_draft.md (Complete manuscript draft)\n")
  cat("   - figure_selection.md (Recommended plots for figures)\n")
  cat("   - summary_statistics.md (Key statistics and results)\n")
  
  return(manuscript_dir)
}

#' Check if RNAseqEr analysis results are available for manuscript generation
#' @description
#' Validates that the necessary files and directories exist for manuscript generation.
#'
#' @param save_dir Directory where analysis results are saved
#' @param dir_lab Label used for the analysis
#'
#' @return List with validation results
#' @export
#'
#' @examples
#' validation <- check_manuscript_prerequisites("analysis_results", "all_celltypes")
#' if(validation$ready) {
#'   generate_manuscript_standalone("analysis_results")
#' }
#'
check_manuscript_prerequisites <- function(save_dir, dir_lab = "all_celltypes") {
  
  # Required directories and files
  required_paths <- list(
    condition_markers = paste0(save_dir, "/outs/", dir_lab, "/tables/condition_mark"),
    cluster_markers = paste0(save_dir, "/outs/", dir_lab, "/tables/cluster_marker"),
    plots = paste0(save_dir, "/outs/", dir_lab, "/plots"),
    reports = paste0(save_dir, "/reports")
  )
  
  # Check existence
  exists <- sapply(required_paths, dir.exists)
  
  # Check for specific files
  files_found <- list()
  
  # Check for condition marker files
  if(exists["condition_markers"]) {
    cond_files <- list.files(required_paths$condition_markers, 
                            pattern = "\\.csv$", 
                            recursive = TRUE, 
                            full.names = TRUE)
    files_found$condition_markers <- length(cond_files) > 0
  } else {
    files_found$condition_markers <- FALSE
  }
  
  # Check for cluster marker files
  if(exists["cluster_markers"]) {
    cluster_files <- list.files(required_paths$cluster_markers, 
                               pattern = "\\.csv$", 
                               recursive = TRUE, 
                               full.names = TRUE)
    files_found$cluster_markers <- length(cluster_files) > 0
  } else {
    files_found$cluster_markers <- FALSE
  }
  
  # Check for plot files
  if(exists["plots"]) {
    plot_files <- list.files(required_paths$plots, 
                            pattern = "\\.pdf$|\\.png$", 
                            recursive = TRUE, 
                            full.names = TRUE)
    files_found$plots <- length(plot_files) > 0
  } else {
    files_found$plots <- FALSE
  }
  
  # Overall readiness
  ready <- all(exists) && all(unlist(files_found))
  
  # Create recommendations
  recommendations <- character()
  if(!exists["condition_markers"]) {
    recommendations <- c(recommendations, "Run find_cond_markers() to generate condition-specific gene lists")
  }
  if(!exists["cluster_markers"]) {
    recommendations <- c(recommendations, "Run int_res_all_mark() or pairwise_dge() to generate cluster markers")
  }
  if(!exists["plots"]) {
    recommendations <- c(recommendations, "Generate plots using plot_list(), save_vln(), or save_feat_plots()")
  }
  if(!exists["reports"]) {
    recommendations <- c(recommendations, "Run generate_report() to create analysis reports")
  }
  
  return(list(
    ready = ready,
    directories_exist = exists,
    files_found = files_found,
    recommendations = recommendations,
    required_paths = required_paths
  ))
}

#' Quick manuscript generation with minimal parameters
#' @description
#' Simplified manuscript generation function that uses sensible defaults
#' and automatically detects available data.
#'
#' @param save_dir Directory where analysis results are saved
#' @param project_title Title for the manuscript
#' @param researcher_name Name of the researcher
#' @param journal_target Target journal (default: "bioRxiv")
#' @param llm_provider LLM provider (default: "openai")
#' @param api_key API key (if NULL, will try to load from environment)
#'
#' @return Path to the manuscript directory
#' @export
#'
#' @examples
#' \dontrun{
#' # Quick manuscript generation
#' manuscript_dir <- quick_manuscript(
#'   save_dir = "analysis_results",
#'   project_title = "My Single-cell Analysis",
#'   researcher_name = "Dr. Smith"
#' )
#' }
#'
quick_manuscript <- function(save_dir,
                           project_title,
                           researcher_name,
                           journal_target = "bioRxiv",
                           llm_provider = "openai",
                           api_key = NULL) {
  
  # Check prerequisites
  validation <- check_manuscript_prerequisites(save_dir)
  
  if(!validation$ready) {
    cat("⚠️  Some prerequisites are missing. Recommendations:\n")
    for(rec in validation$recommendations) {
      cat("   - ", rec, "\n")
    }
    cat("\nYou can still proceed, but the manuscript may be incomplete.\n\n")
  }
  
  # Try to detect conditions of interest from available files
  conditions_of_interest <- NULL
  if(validation$files_found$condition_markers) {
    cond_dir <- paste0(save_dir, "/outs/all_celltypes/tables/condition_mark")
    if(dir.exists(cond_dir)) {
      cond_files <- list.files(cond_dir, recursive = TRUE, pattern = "\\.csv$")
      # Extract condition names from file names
      conditions <- unique(sapply(strsplit(cond_files, "/"), function(x) x[1]))
      conditions_of_interest <- conditions[!is.na(conditions)]
    }
  }
  
  # Generate manuscript with detected parameters
  manuscript_dir <- generate_manuscript_standalone(
    save_dir = save_dir,
    dir_lab = "all_celltypes",
    conditions_of_interest = conditions_of_interest,
    project_title = project_title,
    researcher_name = researcher_name,
    journal_target = journal_target,
    llm_provider = llm_provider,
    api_key = api_key,
    auto_load_env = TRUE
  )
  
  return(manuscript_dir)
} 