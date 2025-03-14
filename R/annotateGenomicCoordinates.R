#' Annotate Multiple Genomic Coordinates Using GTF
#'
#' @description Improved function to annotate genomic coordinates based on a GTF file
#'
#' @param coordinates Data frame with 'chr' and 'pos' columns
#' @param genome Character. Genome build name (e.g., "hg38", "hg19", "mm10")
#' @param gtffile Character. Path to the GTF annotation file
#' @param tss_upstream Numeric. Bases upstream of TSS to define as promoter (default: 2000)
#' @param tss_downstream Numeric. Bases downstream of TSS to include in promoter (default: 200)
#' @param cache_gtf Logical. Whether to cache the GTF data between calls (default: TRUE)
#'
#' @return Data frame with annotation results
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import IRanges
#'
#' @examples
#' # Annotate multiple positions
#' coords <- data.frame(chr = c("chr1", "chr16"), pos = c(1000000, 4334103))
#' results <- annotate_genomic_coordinates(coords, "hg38", "path/to/gencode.v38.annotation.gtf")
# annotate_genomic_coordinates <- function(coordinates, genome, gtffile,
#                                          tss_upstream = 2000, tss_downstream = 200,
#                                          cache_gtf = TRUE) {
  # # Load required packages
  # suppressPackageStartupMessages({
  #   library(rtracklayer)
  #   library(GenomicRanges)
  #   library(IRanges)
  #   library(GenomeInfoDb)
  # })
  #
  # # Validate input
  # if (!is.data.frame(coordinates) || !all(c("chr", "pos") %in% colnames(coordinates))) {
  #   stop("Input must be a data frame with 'chr' and 'pos' columns")
  # }
  #
  # if (!file.exists(gtffile)) {
  #   stop("GTF file not found: ", gtffile)
  # }
  #
  # # Add ID column if not present
  # if (!"id" %in% colnames(coordinates)) {
  #   coordinates$id <- 1:nrow(coordinates)
  # }
  #
  # # Use a cached version of the GTF data if available and requested
  # gtf_env <- new.env()
  # gtf_key <- paste0("gtf_", digest::digest(gtffile))
  #
  # if (cache_gtf && exists(gtf_key, envir = gtf_env)) {
  #   message("Using cached GTF data...")
  #   gtf_data <- get(gtf_key, envir = gtf_env)
  #   genes <- get(paste0(gtf_key, "_genes"), envir = gtf_env)
  # } else {
  #   # Import GTF file only once
  #   message("Importing GTF file...")
  #   gtf_data <- import(gtffile)
  #
  #   # Extract genes
  #   genes <- gtf_data[gtf_data$type == "gene"]
  #
  #   # Ensure genes have a symbol
  #   if ("gene_name" %in% names(mcols(genes))) {
  #     genes$symbol <- genes$gene_name
  #   } else if ("gene_id" %in% names(mcols(genes))) {
  #     genes$symbol <- genes$gene_id
  #   } else {
  #     genes$symbol <- rep(NA, length(genes))
  #   }
  #
  #   # Cache the data if requested
  #   if (cache_gtf) {
  #     assign(gtf_key, gtf_data, envir = gtf_env)
  #     assign(paste0(gtf_key, "_genes"), genes, envir = gtf_env)
  #   }
  # }


  annotate_genomic_coordinates <- function(coordinates, genome, gtffile,
                                           tss_upstream = 2000, tss_downstream = 200,
                                           cache_gtf = TRUE) {
    # Load required packages
    suppressPackageStartupMessages({
      library(rtracklayer)
      library(GenomicRanges)
      library(IRanges)
      library(GenomeInfoDb)
    })

    # Validate input
    if (!is.data.frame(coordinates) || !all(c("chr", "pos") %in% colnames(coordinates))) {
      stop("Input must be a data frame with 'chr' and 'pos' columns")
    }

    # Add ID column if not present
    if (!"id" %in% colnames(coordinates)) {
      coordinates$id <- 1:nrow(coordinates)
    }

    # Use our package-level cache
    gtf_result <- load_cached_gtf(gtffile, force = !cache_gtf)
    gtf_data <- gtf_result$gtf_data
    genes <- gtf_result$genes

  # Create result dataframe
  results <- data.frame(
    id = coordinates$id,
    chr = coordinates$chr,
    position = coordinates$pos,
    genome = rep(genome, nrow(coordinates)),
    location_type = rep(NA, nrow(coordinates)),
    feature = rep(NA, nrow(coordinates)),
    gene_symbol = rep(NA, nrow(coordinates)),
    gene_strand = rep(NA, nrow(coordinates)),
    gene_type = rep(NA, nrow(coordinates)),
    upstream_gene = rep(NA, nrow(coordinates)),
    upstream_distance = rep(NA, nrow(coordinates)),
    downstream_gene = rep(NA, nrow(coordinates)),
    downstream_distance = rep(NA, nrow(coordinates)),
    stringsAsFactors = FALSE
  )

  # First, build chromosome name mapping
  chr_mapping <- character()
  unique_chrs <- unique(coordinates$chr)
  for (chr in unique_chrs) {
    if (chr %in% seqlevels(genes)) {
      chr_mapping[chr] <- chr
    } else if (sub("^chr", "", chr) %in% seqlevels(genes)) {
      chr_mapping[chr] <- sub("^chr", "", chr)
    } else if (paste0("chr", chr) %in% seqlevels(genes)) {
      chr_mapping[chr] <- paste0("chr", chr)
    } else {
      # Set a default if not found
      warning("Chromosome ", chr, " not found in GTF")
      chr_mapping[chr] <- NA
    }
  }

  # Organize genes by chromosome for faster lookup
  genes_by_chr <- list()
  valid_chroms <- seqlevels(genes)
  for (chr in valid_chroms) {
    genes_by_chr[[chr]] <- genes[seqnames(genes) == chr]
  }

  # Process each coordinate
  message("Annotating ", nrow(coordinates), " positions...")

  for (i in 1:nrow(coordinates)) {
    if (i %% 1000 == 0 || i == nrow(coordinates)) {
      message("Processed ", i, " of ", nrow(coordinates), " positions")
    }

    tryCatch({
      chr <- coordinates$chr[i]
      pos <- coordinates$pos[i]

      # Map chromosome name
      chr_to_use <- chr_mapping[chr]
      if (is.na(chr_to_use)) {
        next  # Skip this position if chromosome not found
      }


      # Create query GRange
      queryGR <- GRanges(
        seqnames = chr_to_use,
        ranges = IRanges(start = pos, end = pos)
      )

      # Get genes for this chromosome
      chr_genes <- genes_by_chr[[chr_to_use]]
      if (length(chr_genes) == 0) {
        next  # No genes on this chromosome
      }

      # Check for gene overlaps (is position in a gene?)
      overlaps <- findOverlaps(queryGR, chr_genes, select = "first")

      if (!is.na(overlaps)) {
        # Position is within a gene
        results$location_type[i] <- "genic"
        gene <- chr_genes[overlaps]

        results$gene_symbol[i] <- gene$symbol
        results$gene_strand[i] <- as.character(strand(gene))

        if ("gene_type" %in% names(mcols(gene))) {
          results$gene_type[i] <- gene$gene_type
        } else if ("gene_biotype" %in% names(mcols(gene))) {
          results$gene_type[i] <- gene$gene_biotype
        }

        # Get gene features
        gene_id <- gene$gene_id

        # Get exons, UTRs for this gene (filter directly without re-querying)
        exons <- gtf_data[gtf_data$type == "exon" & gtf_data$gene_id == gene_id]
        five_utrs <- gtf_data[gtf_data$type == "five_prime_utr" & gtf_data$gene_id == gene_id]
        three_utrs <- gtf_data[gtf_data$type == "three_prime_utr" & gtf_data$gene_id == gene_id]

        # Check what feature this position overlaps with
        feature_determined <- FALSE

        if (length(five_utrs) > 0 && !feature_determined) {
          five_utr_overlaps <- findOverlaps(queryGR, five_utrs, select = "first")
          if (!is.na(five_utr_overlaps)) {
            results$feature[i] <- "5' UTR"
            feature_determined <- TRUE
          }
        }

        if (length(three_utrs) > 0 && !feature_determined) {
          three_utr_overlaps <- findOverlaps(queryGR, three_utrs, select = "first")
          if (!is.na(three_utr_overlaps)) {
            results$feature[i] <- "3' UTR"
            feature_determined <- TRUE
          }
        }

        if (length(exons) > 0 && !feature_determined) {
          exon_overlaps <- findOverlaps(queryGR, exons, select = "first")
          if (!is.na(exon_overlaps)) {
            exon <- exons[exon_overlaps]
            if ("exon_number" %in% names(mcols(exon))) {
              results$feature[i] <- paste0("exon ", exon$exon_number)
            } else {
              results$feature[i] <- "exon"
            }
            feature_determined <- TRUE
          }
        }

        # If we get here and feature not determined, position is in an intron
        if (!feature_determined) {
          results$feature[i] <- "intron"
        }
      } else {
        # Position is intergenic
        results$location_type[i] <- "intergenic"

        if (length(chr_genes) > 0) {
          # Find upstream genes (start position is greater than query position)
          upstream_genes <- chr_genes[start(chr_genes) > pos]
          if (length(upstream_genes) > 0) {
            # Find closest upstream gene
            upstream_dists <- start(upstream_genes) - pos
            closest_idx <- which.min(upstream_dists)
            closest_upstream <- upstream_genes[closest_idx]
            results$upstream_gene[i] <- closest_upstream$symbol
            results$upstream_distance[i] <- upstream_dists[closest_idx]
          }

          # Find downstream genes (end position is less than query position)
          downstream_genes <- chr_genes[end(chr_genes) < pos]
          if (length(downstream_genes) > 0) {
            # Find closest downstream gene
            downstream_dists <- pos - end(downstream_genes)
            closest_idx <- which.min(downstream_dists)
            closest_downstream <- downstream_genes[closest_idx]
            results$downstream_gene[i] <- closest_downstream$symbol
            results$downstream_distance[i] <- downstream_dists[closest_idx]
          }
        }
      }
    }, error = function(e) {
      warning("Error annotating position ", i, ": ", e$message)
    })
  }

  message("Annotation complete")
  return(results)
}

# Helper function to annotate a single position
annotate_genomic_position <- function(chr, pos, genome, gtffile,
                                      tss_upstream = 2000, tss_downstream = 200) {
  coords <- data.frame(chr = chr, pos = pos)
  results <- annotate_genomic_coordinates(coords, genome, gtffile,
                                          tss_upstream, tss_downstream)
  return(results[1, ])
}

# Function to pre-load a GTF file for future queries
preload_gtf <- function(gtffile) {
  # Create a digest for the filename as a unique key
  suppressPackageStartupMessages({
    library(digest)
    library(rtracklayer)
    library(GenomicRanges)
  })

  message("Preloading GTF file: ", gtffile)

  # Create environment for caching
  if (!exists("gtf_env", .GlobalEnv)) {
    assign("gtf_env", new.env(), .GlobalEnv)
  }

  gtf_key <- paste0("gtf_", digest(gtffile))

  # Import and process the GTF file
  gtf_data <- import(gtffile)

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

  # Cache the GTF data
  assign(gtf_key, gtf_data, envir = get("gtf_env", .GlobalEnv))
  assign(paste0(gtf_key, "_genes"), genes, envir = get("gtf_env", .GlobalEnv))

  message("GTF file preloaded and ready for use")
  invisible(TRUE)
}
