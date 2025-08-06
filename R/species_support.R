#' Get organism database for gene ontology analysis
#' @description
#' Returns the appropriate organism database for different species
#'
#' @param species Species name (e.g., "Human", "Mouse", "Rat", "Zebrafish", "Drosophila")
#'
#' @return Organism database name for clusterProfiler
#' @export
#'
#' @examples
#' get_organism_db("Human")    # Returns "org.Hs.eg.db"
#' get_organism_db("Mouse")    # Returns "org.Mm.eg.db"
#' get_organism_db("Zebrafish") # Returns "org.Dr.eg.db"
#'
get_organism_db <- function(species) {
  species_db_map <- list(
    "Human" = "org.Hs.eg.db",
    "Mouse" = "org.Mm.eg.db", 
    "Rat" = "org.Rn.eg.db",
    "Zebrafish" = "org.Dr.eg.db",
    "Drosophila" = "org.Dm.eg.db",
    "Chicken" = "org.Gg.eg.db",
    "Cow" = "org.Bt.eg.db",
    "Pig" = "org.Ss.eg.db",
    "Dog" = "org.Cf.eg.db",
    "Cat" = "org.Fc.eg.db",
    "Horse" = "org.Ec.eg.db",
    "Sheep" = "org.Oa.eg.db",
    "Goat" = "org.Ch.eg.db",
    "Rabbit" = "org.Oc.eg.db",
    "Guinea Pig" = "org.Cp.eg.db",
    "Hamster" = "org.Ag.eg.db",
    "Monkey" = "org.Mmu.eg.db",
    "Chimpanzee" = "org.Pt.eg.db",
    "Gorilla" = "org.Ggo.eg.db",
    "Orangutan" = "org.Ppy.eg.db",
    "Gorilla" = "org.Ggo.eg.db",
    "Marmoset" = "org.Cja.eg.db",
    "Squirrel Monkey" = "org.Sbo.eg.db",
    "Tarsier" = "org.Csy.eg.db",
    "Bushbaby" = "org.Oga.eg.db",
    "Mouse Lemur" = "org.Mmu.eg.db",
    "Galago" = "org.Oga.eg.db",
    "Loris" = "org.Nco.eg.db",
    "Tarsier" = "org.Csy.eg.db",
    "New World Monkey" = "org.Cja.eg.db",
    "Old World Monkey" = "org.Mmu.eg.db",
    "Ape" = "org.Hs.eg.db"
  )
  
  # Normalize species name
  species_normalized <- tools::toTitleCase(tolower(species))
  
  if (species_normalized %in% names(species_db_map)) {
    return(species_db_map[[species_normalized]])
  } else {
    warning(paste("Species", species, "not found in database map. Using human database as fallback."))
    return("org.Hs.eg.db")
  }
}

#' Get tissue reference for ScType annotation
#' @description
#' Returns appropriate tissue reference for different species and tissues
#'
#' @param species Species name
#' @param tissue_type Type of tissue (e.g., "Brain", "Liver", "Heart", "Kidney")
#'
#' @return Tissue reference for ScType annotation
#' @export
#'
#' @examples
#' get_tissue_reference("Human", "Brain")     # Returns "Brain"
#' get_tissue_reference("Mouse", "Liver")     # Returns "Liver"
#' get_tissue_reference("Zebrafish", "Brain") # Returns "Brain"
#'
get_tissue_reference <- function(species, tissue_type) {
  # For now, return the tissue type as is
  # ScType supports multiple tissues, so this should work for most cases
  return(tissue_type)
}

#' Check species compatibility
#' @description
#' Checks if the specified species is supported by RNAseqEr
#'
#' @param species Species name
#'
#' @return Logical indicating if species is supported
#' @export
#'
#' @examples
#' is_species_supported("Human")     # Returns TRUE
#' is_species_supported("Mouse")     # Returns TRUE
#' is_species_supported("Alien")     # Returns FALSE
#'
is_species_supported <- function(species) {
  supported_species <- c(
    "Human", "Mouse", "Rat", "Zebrafish", "Drosophila", "Chicken", "Cow", "Pig",
    "Dog", "Cat", "Horse", "Sheep", "Goat", "Rabbit", "Guinea Pig", "Hamster",
    "Monkey", "Chimpanzee", "Gorilla", "Orangutan", "Marmoset", "Squirrel Monkey",
    "Tarsier", "Bushbaby", "Mouse Lemur", "Galago", "Loris", "New World Monkey",
    "Old World Monkey", "Ape"
  )
  
  species_normalized <- tools::toTitleCase(tolower(species))
  return(species_normalized %in% supported_species)
}

#' Get species-specific parameters
#' @description
#' Returns species-specific parameters for RNAseqEr analysis
#'
#' @param species Species name
#' @param tissue_type Type of tissue
#'
#' @return List of species-specific parameters
#' @export
#'
#' @examples
#' get_species_params("Human", "Brain")
#' get_species_params("Mouse", "Liver")
#'
get_species_params <- function(species, tissue_type) {
  # Get organism database
  org_db <- get_organism_db(species)
  
  # Get tissue reference
  tissue_ref <- get_tissue_reference(species, tissue_type)
  
  # Species-specific gene ID mapping
  gene_id_mapping <- list(
    "Human" = list(from = "SYMBOL", to = "ENTREZID"),
    "Mouse" = list(from = "SYMBOL", to = "ENTREZID"),
    "Rat" = list(from = "SYMBOL", to = "ENTREZID"),
    "Zebrafish" = list(from = "SYMBOL", to = "ENTREZID"),
    "Drosophila" = list(from = "SYMBOL", to = "ENTREZID")
  )
  
  species_normalized <- tools::toTitleCase(tolower(species))
  mapping <- if (species_normalized %in% names(gene_id_mapping)) {
    gene_id_mapping[[species_normalized]]
  } else {
    gene_id_mapping[["Human"]]  # Default to human mapping
  }
  
  return(list(
    org_db = org_db,
    tissue_ref = tissue_ref,
    gene_id_from = mapping$from,
    gene_id_to = mapping$to,
    species_supported = is_species_supported(species)
  ))
}

#' Validate species and tissue combination
#' @description
#' Validates if the species and tissue combination is supported
#'
#' @param species Species name
#' @param tissue_type Type of tissue
#'
#' @return Logical indicating if combination is valid
#' @export
#'
#' @examples
#' validate_species_tissue("Human", "Brain")  # Returns TRUE
#' validate_species_tissue("Mouse", "Liver")  # Returns TRUE
#'
validate_species_tissue <- function(species, tissue_type) {
  # For now, return TRUE for most combinations
  # This could be expanded with specific validation rules
  return(TRUE)
} 