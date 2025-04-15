# Define a package-level environment to store loaded data
.pkg_env <- new.env(parent = emptyenv())

#' Get NeoSplice fusion database
#'
#' Retrieves the NeoSplice fusion database, loading it only once when first needed.
#' @keywords internal
#' @return A data frame containing the NeoSplice fusion database
#' @export
get_neosplice_db <- function() {
  if (!exists("neosplice_db", envir = .pkg_env)) {
    message("Loading NeoSplice database for the first time...")
    db_path <- system.file("extdata", "fusions", "NeoSplice_hg38_inframe_fusion.txt.gz",
                           package = "yourpackage")  # Replace "yourpackage" with your actual package name
    .pkg_env$neosplice_db <- parse_neosplice_database(db_path)
  }
  return(.pkg_env$neosplice_db)
}


#' Load and cache GTF data
#'
#' Loads GTF data from a file and caches it for future use.
#' The data is only loaded once per R session.
#'
#' @param gtffile Path to the GTF annotation file
#' @param force Logical. Force reload even if cached (default: FALSE)
#' @return A list containing GTF data and gene information
#' @keywords internal
#' @export
#'
load_cached_gtf <- function(gtffile, force = FALSE) {
  # Check if file exists
  if (!file.exists(gtffile)) {
    stop("GTF file not found: ", gtffile)
  }

  # Create a key for the GTF file
  gtf_key <- paste0("gtf_", digest::digest(gtffile))

  # Check if already cached
  if (!force && exists(gtf_key, envir = .fusion_pkg_env)) {
    message("Using cached GTF data...")
    gtf_data <- get(gtf_key, envir = .fusion_pkg_env)
    genes <- get(paste0(gtf_key, "_genes"), envir = .fusion_pkg_env)
    return(list(gtf_data = gtf_data, genes = genes))
  }

  # Load required packages
  suppressPackageStartupMessages({
    require(rtracklayer)
    require(GenomicRanges)
  })

  # Import GTF file
  message("Importing GTF file...")
  gtf_data <- rtracklayer::import(gtffile)

  # Extract genes
  genes <- gtf_data[gtf_data$type == "gene"]

  # Ensure genes have a symbol
  if ("gene_name" %in% names(mcols(genes))) {
    genes$symbol <- genes$gene_name
  } else if ("gene_id" %in% names(mcols(genes))) {
    genes$symbol <- genes$gene_id
  } else {
    genes$symbol <- rep(NA, length(genes))
  }

  # Cache the data
  assign(gtf_key, gtf_data, envir = .fusion_pkg_env)
  assign(paste0(gtf_key, "_genes"), genes, envir = .fusion_pkg_env)

  return(list(gtf_data = gtf_data, genes = genes))
}
