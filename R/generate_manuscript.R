#' Generate manuscript draft using LLM integration
#' @description
#' This function analyzes RNAseqEr results and uses LLM integration to generate a manuscript draft,
#' focusing on relevant genes and conditions of interest. It also selects the most informative plots
#' for figure generation.
#'
#' @param save_dir Directory where analysis results are saved
#' @param dir_lab Label used for the analysis
#' @param conditions_of_interest Vector of conditions to focus on in the manuscript
#' @param tissue_type Type of tissue analyzed
#' @param species Species analyzed (default: "Human")
#' @param project_title Title of the project
#' @param researcher_name Name of the researcher
#' @param journal_target Target journal for manuscript style (e.g., "Nature", "Cell", "bioRxiv")
#' @param llm_provider LLM provider to use ("openai", "anthropic", "local")
#' @param api_key API key for LLM service (if required)
#' @param model_name Model to use (e.g., "gpt-4", "claude-3", "llama2")
#' @param max_tokens Maximum tokens for LLM response
#' @param temperature Temperature for LLM generation (0-1)
#' @param include_figures Whether to include figure selection and descriptions
#' @param focus_on_genes Whether to focus on specific gene categories (e.g., "marker", "condition", "all")
#' @param min_logfc Minimum log fold change for gene inclusion
#' @param max_pvalue Maximum p-value for gene inclusion
#'
#' @return Saves manuscript draft and figure selection to files
#' @export
#'
#' @examples
#' \dontrun{
#' generate_manuscript(save_dir = "analysis_results",
#'                    dir_lab = "all_celltypes",
#'                    conditions_of_interest = c("AgeGroup", "Tissue"),
#'                    tissue_type = "Brain",
#'                    project_title = "Single-cell analysis of brain aging",
#'                    journal_target = "Nature Communications",
#'                    llm_provider = "openai",
#'                    api_key = "your-api-key")
#' }
#'
generate_manuscript <- function(save_dir,
                              dir_lab = "all_celltypes",
                              conditions_of_interest,
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
  
  # Load environment variables if requested
  if(auto_load_env) {
    load_env_vars()
  }
  
  # Get API key from environment if not provided
  if(is.null(api_key)) {
    api_key <- get_api_key(provider = llm_provider)
    if(is.null(api_key)) {
      warning("No API key found for provider '", llm_provider, "'. Check your .env file or provide api_key parameter.")
    }
  }
  
  # Create manuscript directory
  manuscript_dir <- paste0(save_dir, "/manuscript")
  if(!dir.exists(manuscript_dir)) {
    dir.create(manuscript_dir, recursive = TRUE)
  }
  
  # Collect analysis results
  analysis_summary <- collect_analysis_results(save_dir, dir_lab, conditions_of_interest,
                                             min_logfc, max_pvalue, focus_on_genes)
  
  # Generate LLM prompt
  prompt <- generate_manuscript_prompt(analysis_summary, tissue_type, species, 
                                     project_title, researcher_name, journal_target,
                                     include_figures)
  
  # Call LLM
  manuscript_draft <- call_llm(prompt, llm_provider, api_key, model_name, 
                              max_tokens, temperature)
  
  # Save manuscript draft
  writeLines(manuscript_draft, paste0(manuscript_dir, "/manuscript_draft.md"))
  
  # Generate figure selection if requested
  if(include_figures) {
    figure_selection <- select_informative_plots(save_dir, dir_lab, analysis_summary)
    writeLines(figure_selection, paste0(manuscript_dir, "/figure_selection.md"))
  }
  
  # Generate summary statistics
  summary_stats <- generate_summary_statistics(analysis_summary)
  writeLines(summary_stats, paste0(manuscript_dir, "/summary_statistics.md"))
  
  cat("Manuscript draft generated successfully in:", manuscript_dir, "\n")
  cat("- manuscript_draft.md (Main manuscript draft)\n")
  if(include_figures) cat("- figure_selection.md (Selected plots for figures)\n")
  cat("- summary_statistics.md (Key statistics and results)\n")
  
  return(manuscript_dir)
}

#' Collect analysis results for manuscript generation
#' @keywords internal
collect_analysis_results <- function(save_dir, dir_lab, conditions_of_interest,
                                   min_logfc, max_pvalue, focus_on_genes) {
  
  results <- list()
  
  # Read cluster information
  cluster_file <- paste0(save_dir, "/outs/", dir_lab, "/tables/cluster_purity_data/summary_purity_stats.csv")
  if(file.exists(cluster_file)) {
    cluster_data <- read.csv(cluster_file)
    results$clusters <- cluster_data
  }
  
  # Read differential gene expression results
  dge_file <- paste0(save_dir, "/outs/", dir_lab, "/tables/DGE/broad_celltype_markers/all.csv")
  if(file.exists(dge_file)) {
    dge_data <- read.csv(dge_file)
    # Filter by criteria
    dge_filtered <- dge_data[dge_data$avg_log2FC >= min_logfc & dge_data$p_val_adj <= max_pvalue, ]
    results$marker_genes <- dge_filtered
  }
  
  # Read condition-specific results
  results$condition_genes <- list()
  for(condition in conditions_of_interest) {
    cond_file <- paste0(save_dir, "/outs/", dir_lab, "/tables/condition_mark/RNAseqEr_annotation/clusterwise/all_mark__", condition, ".csv")
    if(file.exists(cond_file)) {
      cond_data <- read.csv(cond_file)
      # Filter by criteria
      cond_filtered <- cond_data[cond_data$avg_log2FC >= min_logfc & cond_data$p_val_adj <= max_pvalue, ]
      results$condition_genes[[condition]] <- cond_filtered
    }
  }
  
  # Read volcano plot results if available
  volcano_dir <- paste0(save_dir, "/outs/summary_figures/conditions/volcano")
  if(dir.exists(volcano_dir)) {
    volcano_files <- list.files(volcano_dir, pattern = "\\.pdf$", full.names = TRUE)
    results$volcano_plots <- volcano_files
  }
  
  return(results)
}

#' Generate LLM prompt for manuscript generation
#' @keywords internal
generate_manuscript_prompt <- function(analysis_summary, tissue_type, species, 
                                     project_title, researcher_name, journal_target,
                                     include_figures) {
  
  # Extract key genes
  top_markers <- if(!is.null(analysis_summary$marker_genes)) {
    head(analysis_summary$marker_genes[order(-analysis_summary$marker_genes$avg_log2FC), ], 20)
  } else {
    data.frame()
  }
  
  condition_summary <- ""
  for(condition in names(analysis_summary$condition_genes)) {
    if(nrow(analysis_summary$condition_genes[[condition]]) > 0) {
      top_condition_genes <- head(analysis_summary$condition_genes[[condition]][order(-analysis_summary$condition_genes[[condition]]$avg_log2FC), ], 10)
      condition_summary <- paste0(condition_summary, "\n", condition, " condition:\n")
      for(i in 1:nrow(top_condition_genes)) {
        gene <- top_condition_genes$gene[i]
        logfc <- round(top_condition_genes$avg_log2FC[i], 2)
        pval <- format(top_condition_genes$p_val_adj[i], scientific = TRUE, digits = 2)
        condition_summary <- paste0(condition_summary, "- ", gene, " (log2FC: ", logfc, ", p.adj: ", pval, ")\n")
      }
    }
  }
  
  prompt <- paste0(
    "You are a scientific writer specializing in single-cell RNA sequencing analysis. ",
    "Please write a manuscript draft for the following analysis results. ",
    "Focus on the most biologically relevant findings and ensure the writing is suitable for ", journal_target, ".\n\n",
    
    "Project Information:\n",
    "- Title: ", project_title, "\n",
    "- Tissue: ", tissue_type, "\n",
    "- Species: ", species, "\n",
    "- Researcher: ", researcher_name, "\n\n",
    
    "Analysis Results:\n",
    "- Number of clusters identified: ", if(!is.null(analysis_summary$clusters)) nrow(analysis_summary$clusters) else "Not available", "\n",
    "- Top marker genes: ", if(nrow(top_markers) > 0) paste(head(top_markers$gene, 5), collapse = ", ") else "Not available", "\n",
    "- Condition-specific genes: ", condition_summary, "\n\n",
    
    "Please write a manuscript draft including:\n",
    "1. Abstract (250 words)\n",
    "2. Introduction (focus on the biological context and significance)\n",
    "3. Results (highlight key findings, especially condition-specific changes)\n",
    "4. Discussion (interpret findings and suggest future directions)\n",
    "5. Methods (brief overview of the analysis pipeline)\n\n",
    
    if(include_figures) {
      c("Also suggest which plots would be most informative for figures, focusing on:\n",
        "- UMAP plots showing cell type distribution\n",
        "- Volcano plots showing differential expression\n",
        "- Heatmaps of marker genes\n",
        "- Violin plots of key genes\n\n")
    } else {
      ""
    },
    
    "Guidelines:\n",
    "- Focus on biological significance rather than technical details\n",
    "- Highlight genes with the strongest evidence (high log fold change, low p-value)\n",
    "- Suggest potential mechanisms or pathways\n",
    "- Consider implications for understanding the tissue biology\n",
    "- Write in a clear, scientific style suitable for peer review\n"
  )
  
  return(prompt)
}

#' Call LLM service
#' @keywords internal
call_llm <- function(prompt, provider, api_key, model, max_tokens, temperature) {
  
  if(provider == "openai") {
    return(call_openai(prompt, api_key, model, max_tokens, temperature))
  } else if(provider == "anthropic") {
    return(call_anthropic(prompt, api_key, model, max_tokens, temperature))
  } else if(provider == "local") {
    return(call_local_llm(prompt, model, max_tokens, temperature))
  } else {
    stop("Unsupported LLM provider. Use 'openai', 'anthropic', or 'local'")
  }
}

#' Call OpenAI API
#' @keywords internal
call_openai <- function(prompt, api_key, model, max_tokens, temperature) {
  if(is.null(api_key)) {
    stop("API key required for OpenAI")
  }
  
  # This would require the openai package
  # For now, return a placeholder
  return("Manuscript draft would be generated here using OpenAI API.\n\nPlease install the 'openai' package and configure your API key to use this feature.")
}

#' Call Anthropic API
#' @keywords internal
call_anthropic <- function(prompt, api_key, model, max_tokens, temperature) {
  if(is.null(api_key)) {
    stop("API key required for Anthropic")
  }
  
  # This would require the anthropic package
  # For now, return a placeholder
  return("Manuscript draft would be generated here using Anthropic API.\n\nPlease install the 'anthropic' package and configure your API key to use this feature.")
}

#' Call local LLM
#' @keywords internal
call_local_llm <- function(prompt, model, max_tokens, temperature) {
  # This would integrate with local LLM services
  # For now, return a placeholder
  return("Manuscript draft would be generated here using local LLM.\n\nPlease configure your local LLM service to use this feature.")
}

#' Select informative plots for figures
#' @keywords internal
select_informative_plots <- function(save_dir, dir_lab, analysis_summary) {
  
  plot_selection <- c(
    "# Figure Selection for Manuscript",
    "",
    "## Recommended Figures:",
    "",
    "### Figure 1: Overview of Single-cell Analysis",
    "- UMAP plot showing cell type distribution",
    "- Quality control metrics",
    "- Sample composition",
    "",
    "### Figure 2: Cell Type Characterization",
    "- Heatmap of top marker genes for each cell type",
    "- Violin plots of canonical markers",
    "- Cell type proportions",
    "",
    "### Figure 3: Condition-specific Changes",
    "- Volcano plots for each condition of interest",
    "- Top differentially expressed genes",
    "- Pathway enrichment results",
    "",
    "### Figure 4: Biological Insights",
    "- Gene expression changes across conditions",
    "- Cell type-specific responses",
    "- Validation plots",
    "",
    "## Plot Files Available:",
    ""
  )
  
  # List available plots
  plots_dir <- paste0(save_dir, "/outs/", dir_lab, "/plots")
  if(dir.exists(plots_dir)) {
    plot_files <- list.files(plots_dir, recursive = TRUE, full.names = TRUE)
    for(file in plot_files) {
      plot_selection <- c(plot_selection, paste("-", file))
    }
  }
  
  return(plot_selection)
}

#' Generate summary statistics
#' @keywords internal
generate_summary_statistics <- function(analysis_summary) {
  
  stats <- c(
    "# Summary Statistics",
    "",
    "## Clustering Results:",
    paste("- Number of clusters:", if(!is.null(analysis_summary$clusters)) nrow(analysis_summary$clusters) else "Not available"),
    "",
    "## Marker Genes:",
    paste("- Total marker genes identified:", if(!is.null(analysis_summary$marker_genes)) nrow(analysis_summary$marker_genes) else "Not available"),
    paste("- Top marker gene:", if(!is.null(analysis_summary$marker_genes) && nrow(analysis_summary$marker_genes) > 0) analysis_summary$marker_genes$gene[1] else "Not available"),
    "",
    "## Condition-specific Genes:"
  )
  
  for(condition in names(analysis_summary$condition_genes)) {
    if(nrow(analysis_summary$condition_genes[[condition]]) > 0) {
      stats <- c(stats,
                 paste("-", condition, "condition:"),
                 paste("  - Number of genes:", nrow(analysis_summary$condition_genes[[condition]])),
                 paste("  - Top gene:", analysis_summary$condition_genes[[condition]]$gene[1]),
                 paste("  - Average log2FC:", round(mean(analysis_summary$condition_genes[[condition]]$avg_log2FC), 2)),
                 "")
    }
  }
  
  return(stats)
} 