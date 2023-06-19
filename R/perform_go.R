#' perfome_go
#' @description
#' Run gene ontology using cluster profiler
#'
#' @param srt Seurat object
#' @param gene_list list of genes to use for gene ontology analyses
#' @param min_log2FC log2 fold cutoff for genes to be considered in gene ontology
#' analysis
#' @param reverse if genes smaller than the log2 FC cuttoff should be considered
#' instead, for example for smaler than -0.25
#' @param translate_gene_id_from current format of gene ids, default is "SYMBOL"
#' @param translate_gene_id_to specify which gene id's should be used (normally
#' ENTREZID), default is c("ENTREZID", "PATH", "GO", "ALIAS", "GENENAME")
#' @param org_use specify reference for your organism, default is human ("org.Dr.eg.db")
#' @param ontology what kind of gene ontology should be performed. Default is "BP" for
#' biological pathway but "MF" for molecular function for example is valid, too.
#' @param pvalue_cutoff p value cutoff for significant gene ontology results
#' @param qvalue_cutoff q value cutoff for significant gene ontology results
#' @param read_able gene names in GO output readable? Default is TRUE
#' @import clusterProfiler
#'
#' @return returns a gene ontology result (add data format) ready for plotting
#' and saving
#' @export
#'
#' @examples
#' library(Seurat)
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' data(dge_1_0)
#'
#' perform_go(pbmc_small, gene_list = dge_1_0)
#'
perform_go <- function(srt,
                       gene_list,
                       min_log2FC = 0.25,
                       reverse = FALSE,
                       translate_gene_id_from = "SYMBOL",
                       translate_gene_id_to = "ENTREZID",
                       org_use = "org.Hs.eg.db",
                       ontology = "BP",
                       pvalue_cutoff = 0.05,
                       qvalue_cutoff = 0.05,
                       read_able = TRUE) {
  if (reverse == TRUE) {
    genes_inc <- as.data.frame(subset(gene_list, gene_list$avg_log2FC < min_log2FC))
    # allows testeing for genes more expressed in the other group
  } else {
    genes_inc <- as.data.frame(subset(gene_list, gene_list$avg_log2FC > min_log2FC))
  }

  if (class(genes_inc$X) != "NULL") { # this happens when reading in a saves .csv
    genes <- genes_inc$X
  } else if (class(genes_inc$gene) != "NULL") {
    genes <- genes_inc$gene
  } else {
    genes <- rownames(gene_list) # this is for newly generated gene lists
  }
  genes_entrez <- bitr(genes,
    fromType = translate_gene_id_from,
    toType = translate_gene_id_to,
    OrgDb = org_use
  )
  go_genes <- genes_entrez$ENTREZID
  ego <- enrichGO(
    gene = go_genes,
    universe = names(rownames(srt)),
    OrgDb = org_use,
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = pvalue_cutoff,
    qvalueCutoff = qvalue_cutoff,
    readable = read_able
  )
  return(ego)
}
