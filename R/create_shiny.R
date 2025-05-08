#' Make Shiny app using ShinyCell
#'
#' @param seur_obj seurat object
#' @param shiny_name name used in Shiny tab if only one seurat object is provided.
#' Default is "all_celltypes".
#' @param read_file TRUE/ FALSE whether Seurat objects are read from provided file
#' paths which are up to 2 locations provided in file_dir_1 and file_dir_2.
#' @param file_dir_1 file path to first seurat object saved as RDS file
#' that contains highest level information as in all celltypes together.
#' Needs to be provided if read_file is set to TRUE. Default is FALSE.
#' @param file_dir_2 file path to cell lineage seurat objects. Default is FALSE.
#' @param ext_pattern Extension pattern of saved seurat object files. Default is
#' ".RDS" while ".rds" works as well.
#' @param default_1 First metadata column that should be displayed when opening shiny
#' app.
#' @param default_2 Second metadata column that should be displayed when opening shiny
#' app.
#' @param assay_use Seurat objects: "RNA" or "integrated" assay, default is "RNA"
#' @param gex_slot slot in single-cell assay to plot. Default is to use the "data" slot
#' @param gene_mapping specifies whether to convert human / mouse Ensembl gene IDs
#' (e.g. ENSG000xxx / ENSMUSG000xxx) into "user-friendly" gene symbols.
#' Set this to TRUE if you are using Ensembl gene IDs.
#' Default is FALSE which is not to perform any conversion.
#' Alternatively, users can supply a named vector where names(gene.mapping)
#' correspond to the actual gene identifiers in the gene expression matrix
#' and gene.mapping correspond to new identifiers to map to.
#' @param default_gene1 provide first default gene to be shown in Shiny app.
#' @param default_gene2 provide second default gene to be shown in Shiny app.
#' @param save_dir directory in which to save Shiny output files.
#' @param default_multigene character vector of genes that should be shown as
#' default in bubbleplot and heatmap
#' @param default_dimred "UMAP" (default) or "TSNE" as embeddings for plots.
#' @param author Primary author name of study. Default is "TBC".
#' @param title Title of study. Default is "TBC".
#' @param journal journal in which study is published. Default is "TBC".
#' @param volume Citation information, journal volume. Default is "TBC".
#' @param page Citation information, journal page number Default is "TBC".
#' @param year Citation information, Publication year. Default is "TBC".
#' @param doi Citation information, DOI. Default is "TBC".
#' @param link Link to manuscript/paper. Default is "TBC".
#' @param shiny_title title of shiny app.
#'
#' @return server.R and ui.R required for shiny app
#' @import ShinyCell shinyhelper
#' @export
#'
#' @examples
#' create_shiny(cns)
create_shiny <- function(seur_obj,
                        shiny_name = "all_celltypes",
                        read_file = FALSE,
                        file_dir_1 = FALSE,
                        file_dir_2 = FALSE,
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
                        ){
  prefix_list <- list()

  if(read_file == FALSE){
  scConf1 = createConfig(seur_obj,
                         meta.to.include = colnames(seur_obj@meta.data))

  scConf1 = modDefault(scConf1, default1 = default_1 , default2 = default_2)

  shiny_dir <- paste0(save_dir, "shiny_app/")

  if(dir.exists(shiny_dir) == FALSE){
    dir.create(shiny_dir)
  }

  makeShinyFiles(seur_obj,
                 scConf1,
                 gex.assay = assay_use,
                 gex.slot = gex_slot,
                 gene.mapping = gene_mapping,
                 shiny.prefix = shiny_name,
                 shiny.dir = shiny_dir ,
                 default.gene1 = default_gene1,
                 default.gene2 = default_gene2,
                 default.multigene = default_multigene,
                 default.dimred = default_dimred)

  file_names <- shiny_name
  
  prefix_list <- shiny_name
  }

  if(read_file == TRUE){
    if(file_dir_1 != FALSE){
    file_names <- list.files(file_dir_1, pattern = ext_pattern)
    file_path <- paste(file_dir_1, file_names, sep = "/")
    }else{
      print("Please provide file path to seurat objects.")
    }

    if(file_dir_2 != FALSE){
    file_names_2 <- list.files(file_dir_2, pattern = ext_pattern)
    file_names <- c(file_names, file_names_2)
    file_path_2 <- paste(file_dir_2, file_names_2, sep = "/")
    file_path = c(file_path, file_path_2)
    }

    for(i in 1: length(file_names)){
      curr_seur <- readRDS(file_path[i])
      curr_name <- file_names[i]

      print(curr_seur)

      curr_name_save <- strsplit(curr_name, ".R")[[1]][1]


      scConf1 = createConfig(curr_seur,
                             meta.to.include = colnames(curr_seur@meta.data))

      scConf1 = modDefault(scConf1, default1 = default_1 , default2 = default_2)

      shiny_dir <- paste0(save_dir, "/shiny_app/")

      if(dir.exists(shiny_dir) == FALSE){
        dir.create(shiny_dir)
      }

      makeShinyFiles(curr_seur,
                     scConf1,
                     gex.assay = assay_use,
                     gex.slot = gex_slot,
                     gene.mapping = gene_mapping,
                     shiny.prefix = curr_name_save,
                     shiny.dir = shiny_dir ,
                     default.gene1 = default_gene1,
                     default.gene2 = default_gene2,
                     default.multigene = default_multigene,
                     default.dimred = default_dimred)

      prefix_list[i] <- curr_name_save

    }

  }

  citation = list(
    author  = author,
    title   = title,
    journal = journal,
    volume  = volume,
    page    = page,
    year    = year,
    doi     = doi,
    link    = link)

  makeShinyCodesMulti(
    shiny.title = shiny_title,
    shiny.footnotes = citation,
    shiny.prefix = prefix_list,
    shiny.headers = prefix_list,
    shiny.dir = shiny_dir)


}













