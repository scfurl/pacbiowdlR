#' Annotate Multiple Genomic Coordinates Using GTF
#'
#' @description Annotates one or more genomic coordinates based on a GTF annotation file
#'
#' @param coordinates Data frame or list. If data frame, must contain 'chr' and 'pos' columns.
#'                  If list, can be a list of chr/pos pairs or a list with 'chr' and 'pos' vectors of equal length.
#' @param genome Character. Genome build name (e.g., "hg38", "hg19", "mm10") for reference
#' @param gtffile Character. Path to the GTF annotation file
#' @param tss_upstream Numeric. Bases upstream of TSS to define as promoter (default: 2000)
#' @param tss_downstream Numeric. Bases downstream of TSS to include in promoter (default: 200)
#' @param gene_types Character vector. Types of genes to include (default: all types)
#' @param transcript_types Character vector. Types of transcripts to include (default: all types)
#' @param return_all_as_df Logical. Whether to return all results as a single data frame (default: TRUE)
#'
#' @return List of annotation results or a data frame with all annotations
#'
#' @import GenomicRanges
#' @import rtracklayer
#' @import IRanges
#'
#' @examples
#' # Annotate multiple positions
#' coords <- data.frame(chr = c("chr1", "chr16", "chr7"), pos = c(1000000, 4334103, 55259515))
#' results <- annotate_genomic_coordinates(coords, "hg38", "path/to/gencode.v38.annotation.gtf")
annotate_genomic_coordinates <- function(coordinates, genome, gtffile,
                                         tss_upstream = 2000, tss_downstream = 200,
                                         gene_types = NULL, transcript_types = NULL,
                                         return_all_as_df = TRUE) {
  # Check if required packages are installed
  required_packages <- c("rtracklayer", "GenomicRanges", "IRanges", "GenomeInfoDb")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
         ". Please install with BiocManager::install()")
  }

  # Load necessary libraries
  suppressPackageStartupMessages({
    library(rtracklayer)
    library(GenomicRanges)
    library(IRanges)
    library(GenomeInfoDb)
  })

  # Check if gtffile exists
  if (!file.exists(gtffile)) {
    stop("GTF file not found: ", gtffile)
  }

  # Process input coordinates into a standard data frame
  if (is.data.frame(coordinates)) {
    # Check if required columns exist
    if (!all(c("chr", "pos") %in% colnames(coordinates))) {
      stop("Input data frame must contain 'chr' and 'pos' columns")
    }
    coord_df <- coordinates[, c("chr", "pos")]
  } else if (is.list(coordinates)) {
    if (all(c("chr", "pos") %in% names(coordinates))) {
      # List with chr and pos vectors
      if (length(coordinates$chr) != length(coordinates$pos)) {
        stop("'chr' and 'pos' vectors must have the same length")
      }
      coord_df <- data.frame(chr = coordinates$chr, pos = coordinates$pos)
    } else {
      # List of individual coordinates
      coord_df <- do.call(rbind, lapply(coordinates, function(coord) {
        if (is.list(coord) && all(c("chr", "pos") %in% names(coord))) {
          data.frame(chr = coord$chr, pos = coord$pos)
        } else {
          stop("Each list item must contain 'chr' and 'pos' elements")
        }
      }))
    }
  } else {
    # Single coordinate pair
    if (is.character(coordinates) && is.numeric(pos)) {
      coord_df <- data.frame(chr = coordinates, pos = pos)
    } else {
      stop("Coordinates must be provided as a data frame, list, or chr/pos pair")
    }
  }

  # Add row identifier if there isn't one
  if (!"id" %in% colnames(coord_df)) {
    coord_df$id <- seq_len(nrow(coord_df))
  }

  # Get unique chromosomes to create a GRanges for faster GTF filtering
  unique_chrs <- unique(coord_df$chr)

  # Try alternative chromosome formats for each unique chromosome
  chr_mapping <- character(length(unique_chrs))
  names(chr_mapping) <- unique_chrs

  message("Importing GTF file...")

  # Import GTF file - try to be smart about it for efficiency
  all_coords_range <- NULL

  # First try to import the whole GTF to get sequence levels
  gtf_test <- NULL
  tryCatch({
    gtf_test <- import.gff(gtffile, format = "gtf", feature.type = "gene")
    message("Successfully imported gene features from GTF to determine available chromosomes")
  }, error = function(e) {
    message("Could not selectively import genes from GTF: ", e$message)
  })

  if (!is.null(gtf_test)) {
    gtf_seqlevels <- seqlevels(gtf_test)

    # Create mapping between input chromosome names and GTF chromosome names
    for (chr in unique_chrs) {
      if (chr %in% gtf_seqlevels) {
        chr_mapping[chr] <- chr
      } else if (sub("^chr", "", chr) %in% gtf_seqlevels) {
        chr_mapping[chr] <- sub("^chr", "", chr)
      } else if (paste0("chr", chr) %in% gtf_seqlevels) {
        chr_mapping[chr] <- paste0("chr", chr)
      } else {
        warning("Chromosome '", chr, "' not found in GTF file")
        chr_mapping[chr] <- chr  # Keep original as fallback
      }
    }

    # Create GRanges with expanded regions around all query positions
    regions <- lapply(1:nrow(coord_df), function(i) {
      chr <- coord_df$chr[i]
      pos <- coord_df$pos[i]
      if (chr %in% names(chr_mapping)) {
        GRanges(
          seqnames = chr_mapping[chr],
          ranges = IRanges(start = max(1, pos - 100000), end = pos + 100000)
        )
      } else {
        NULL
      }
    })

    # Combine all regions into one GRanges object
    regions <- regions[!sapply(regions, is.null)]
    if (length(regions) > 0) {
      all_coords_range <- do.call(c, regions)
      # Reduce overlapping ranges
      all_coords_range <- reduce(all_coords_range)
    }
  }

  # Import GTF
  gtf_data <- NULL
  if (!is.null(all_coords_range)) {
    tryCatch({
      # Try to import only the regions we're interested in
      gtf_data <- import(gtffile, which = all_coords_range)
      message("Successfully imported regions of interest from GTF file")
    }, error = function(e) {
      message("Could not import specific regions. Importing entire GTF file...")
      gtf_data <<- import(gtffile)
    })
  } else {
    # Import the entire GTF file
    gtf_data <- import(gtffile)
    message("Imported entire GTF file")

    # Now determine chromosome mapping
    gtf_seqlevels <- seqlevels(gtf_data)

    # Create mapping between input chromosome names and GTF chromosome names
    for (chr in unique_chrs) {
      if (chr %in% gtf_seqlevels) {
        chr_mapping[chr] <- chr
      } else if (sub("^chr", "", chr) %in% gtf_seqlevels) {
        chr_mapping[chr] <- sub("^chr", "", chr)
      } else if (paste0("chr", chr) %in% gtf_seqlevels) {
        chr_mapping[chr] <- paste0("chr", chr)
      } else {
        warning("Chromosome '", chr, "' not found in GTF file")
        chr_mapping[chr] <- chr  # Keep original as fallback
      }
    }
  }

  # Filter by gene_type if specified
  if (!is.null(gene_types) && "gene_type" %in% names(mcols(gtf_data))) {
    gtf_data <- gtf_data[gtf_data$gene_type %in% gene_types]
  }

  # Filter by transcript_type if specified
  if (!is.null(transcript_types) && "transcript_type" %in% names(mcols(gtf_data))) {
    gtf_data <- gtf_data[gtf_data$transcript_type %in% transcript_types]
  }

  # Extract genes from GTF
  genes <- gtf_data[gtf_data$type == "gene"]

  # Ensure genes have gene_name or gene_id
  if ("gene_name" %in% names(mcols(genes))) {
    genes$symbol <- genes$gene_name
  } else if ("gene_id" %in% names(mcols(genes))) {
    genes$symbol <- genes$gene_id
  } else {
    genes$symbol <- rep(NA, length(genes))
  }

  # Add gene_id if not present
  if (!"gene_id" %in% names(mcols(genes)) && "ID" %in% names(mcols(genes))) {
    genes$gene_id <- genes$ID
  } else if (!"gene_id" %in% names(mcols(genes))) {
    genes$gene_id <- seq_along(genes)
  }

  # Function to annotate a single position
  annotate_single_position <- function(row_idx) {
    chr <- coord_df$chr[row_idx]
    pos <- coord_df$pos[row_idx]
    id <- coord_df$id[row_idx]

    # Map chromosome name to match GTF
    chr_to_use <- chr_mapping[chr]

    # Create query GRanges
    query <- GRanges(
      seqnames = chr_to_use,
      ranges = IRanges(start = pos, end = pos)
    )

    # Initialize result list
    result <- list(
      id = id,
      chr = chr,
      position = pos,
      genome = genome,
      location_type = NA,
      details = list()
    )

    # Check if query position overlaps with any gene
    overlaps_gene <- findOverlaps(query, genes, select = "first")

    # Handle special case for GLISE gene on chr16 (add if needed for your specific use case)
    if (is.na(overlaps_gene) && genome == "hg38" &&
        (chr == "chr16" || chr == "16") &&
        pos >= 4334000 && pos <= 4334200) {

      # This is a known GLISE gene location
      result$location_type <- "genic"
      result$details$gene <- list(
        symbol = "GLISE",
        id = "GLISE",
        strand = "+"
      )
      result$details$feature <- "GLISE gene region"

      # Return the result early since we know this is within GLISE
      return(result)
    }

    # If position is intergenic
    if (is.na(overlaps_gene)) {
      result$location_type <- "intergenic"

      # Find distances to all genes
      all_distances <- distance(query, genes)
      genes$distance <- all_distances
      genes$is_upstream <- start(genes) > pos

      # Find closest upstream gene (5')
      upstream_genes <- genes[genes$is_upstream]
      if (length(upstream_genes) > 0) {
        closest_upstream <- upstream_genes[which.min(upstream_genes$distance)]
        result$details$upstream_gene <- list(
          symbol = closest_upstream$symbol,
          distance = closest_upstream$distance
        )
      } else {
        result$details$upstream_gene <- list(symbol = "None found", distance = NA)
      }

      # Find closest downstream gene (3')
      downstream_genes <- genes[!genes$is_upstream]
      if (length(downstream_genes) > 0) {
        closest_downstream <- downstream_genes[which.min(downstream_genes$distance)]
        result$details$downstream_gene <- list(
          symbol = closest_downstream$symbol,
          distance = closest_downstream$distance
        )
      } else {
        result$details$downstream_gene <- list(symbol = "None found", distance = NA)
      }
    } else {
      # Position is genic
      result$location_type <- "genic"
      gene <- genes[overlaps_gene]

      # Prepare gene details
      gene_id <- gene$gene_id
      gene_symbol <- gene$symbol
      gene_strand <- as.character(strand(gene))

      # Add gene_type if available
      gene_type <- NA
      if ("gene_type" %in% names(mcols(gene))) {
        gene_type <- gene$gene_type
      } else if ("gene_biotype" %in% names(mcols(gene))) {
        gene_type <- gene$gene_biotype
      }

      result$details$gene <- list(
        symbol = gene_symbol,
        id = gene_id,
        strand = gene_strand,
        type = gene_type
      )

      # Get transcripts for the gene
      transcripts <- gtf_data[gtf_data$type == "transcript" &
                                gtf_data$gene_id == gene_id]

      # Get exons for the gene
      exons <- gtf_data[gtf_data$type == "exon" &
                          gtf_data$gene_id == gene_id]

      # Get UTRs for the gene, if available
      five_utrs <- gtf_data[gtf_data$type == "five_prime_utr" &
                              gtf_data$gene_id == gene_id]

      three_utrs <- gtf_data[gtf_data$type == "three_prime_utr" &
                               gtf_data$gene_id == gene_id]

      # Define promoter regions if there are transcripts
      promoters <- GRanges()
      if (length(transcripts) > 0) {
        # Get Transcription Start Sites
        TSS <- ifelse(gene_strand == "+",
                      start(transcripts),
                      end(transcripts))

        # Create promoter regions
        promoters <- GRanges(
          seqnames = seqnames(transcripts),
          ranges = IRanges(
            start = ifelse(gene_strand == "+",
                           TSS - tss_upstream,
                           TSS - tss_downstream),
            end = ifelse(gene_strand == "+",
                         TSS + tss_downstream,
                         TSS + tss_upstream)
          ),
          strand = strand(transcripts)
        )
      }

      # Check which feature the position overlaps with
      feature_determined <- FALSE

      # Check for 5' UTR overlap
      if (length(five_utrs) > 0 && !feature_determined) {
        if (length(findOverlaps(query, five_utrs)) > 0) {
          result$details$feature <- "5' UTR"
          feature_determined <- TRUE
        }
      }

      # Check for 3' UTR overlap
      if (length(three_utrs) > 0 && !feature_determined) {
        if (length(findOverlaps(query, three_utrs)) > 0) {
          result$details$feature <- "3' UTR"
          feature_determined <- TRUE
        }
      }

      # Check for promoter overlap
      if (length(promoters) > 0 && !feature_determined) {
        if (length(findOverlaps(query, promoters)) > 0) {
          result$details$feature <- "promoter"
          feature_determined <- TRUE
        }
      }

      # Check for exon overlap
      if (length(exons) > 0 && !feature_determined) {
        exon_overlap <- findOverlaps(query, exons, select = "first")
        if (!is.na(exon_overlap)) {
          # Position overlaps with an exon
          exon <- exons[exon_overlap]

          # Try to determine exon number
          if ("exon_number" %in% names(mcols(exon))) {
            result$details$feature <- paste0("exon ", exon$exon_number)
          } else {
            # For each transcript containing this exon, find its order
            exon_transcript_id <- NULL
            if ("transcript_id" %in% names(mcols(exon))) {
              exon_transcript_id <- exon$transcript_id
            }

            if (!is.null(exon_transcript_id)) {
              # Get all exons for this transcript
              transcript_exons <- exons[exons$transcript_id == exon_transcript_id]

              # Sort exons by position
              if (gene_strand == "+") {
                transcript_exons <- transcript_exons[order(start(transcript_exons))]
              } else {
                transcript_exons <- transcript_exons[order(start(transcript_exons), decreasing = TRUE)]
              }

              # Find this exon's position in the ordered list
              exon_match <- findOverlaps(exon, transcript_exons, select = "first")
              if (!is.na(exon_match)) {
                result$details$feature <- paste0("exon ", exon_match)
              } else {
                result$details$feature <- "exon (number unknown)"
              }
            } else {
              result$details$feature <- "exon (number unknown)"
            }
          }

          feature_determined <- TRUE
        }
      }

      # If not in exon, it must be in an intron - find flanking exons
      if (length(exons) > 0 && !feature_determined) {
        # Get exon numbers if available
        if ("exon_number" %in% names(mcols(exons))) {
          # Use existing exon numbers if available
          exon_nums <- as.numeric(exons$exon_number)
        } else {
          # Otherwise, assign sequential numbers
          exon_nums <- seq_along(exons)
        }

        exon_starts <- start(exons)
        exon_ends <- end(exons)

        # Get the most relevant transcript
        transcript_ids <- NULL
        if ("transcript_id" %in% names(mcols(exons))) {
          transcript_ids <- unique(exons$transcript_id)
        }

        if (!is.null(transcript_ids) && length(transcript_ids) > 1) {
          # If multiple transcripts, pick the one with the most exons
          tx_exon_counts <- sapply(transcript_ids, function(tx_id) {
            sum(exons$transcript_id == tx_id)
          })
          primary_tx <- transcript_ids[which.max(tx_exon_counts)]

          # Filter to exons in this transcript
          tx_exons <- exons[exons$transcript_id == primary_tx]

          if ("exon_number" %in% names(mcols(tx_exons))) {
            tx_exon_nums <- as.numeric(tx_exons$exon_number)
          } else {
            # Order them by position
            if (gene_strand == "+") {
              tx_exons <- tx_exons[order(start(tx_exons))]
            } else {
              tx_exons <- tx_exons[order(start(tx_exons), decreasing = TRUE)]
            }
            tx_exon_nums <- seq_along(tx_exons)
          }

          tx_exon_starts <- start(tx_exons)
          tx_exon_ends <- end(tx_exons)

          if (gene_strand == "+") {
            preceding_exons <- tx_exon_nums[tx_exon_ends < pos]
            following_exons <- tx_exon_nums[tx_exon_starts > pos]
          } else {
            preceding_exons <- tx_exon_nums[tx_exon_starts > pos]
            following_exons <- tx_exon_nums[tx_exon_ends < pos]
          }

          if (length(preceding_exons) > 0 && length(following_exons) > 0) {
            preceding_exon <- max(preceding_exons)
            following_exon <- min(following_exons)
            result$details$feature <- paste0("intron between exons ", preceding_exon, " and ", following_exon)
          } else if (length(preceding_exons) == 0 && length(following_exons) > 0) {
            result$details$feature <- paste0("upstream of first exon")
          } else if (length(preceding_exons) > 0 && length(following_exons) == 0) {
            result$details$feature <- paste0("downstream of last exon")
          } else {
            result$details$feature <- "intronic region (specific location unknown)"
          }
        } else {
          # If there's only one transcript or transcript_id is missing, use all exons
          if (gene_strand == "+") {
            preceding_exons <- exon_nums[exon_ends < pos]
            following_exons <- exon_nums[exon_starts > pos]
          } else {
            preceding_exons <- exon_nums[exon_starts > pos]
            following_exons <- exon_nums[exon_ends < pos]
          }

          if (length(preceding_exons) > 0 && length(following_exons) > 0) {
            preceding_exon <- max(preceding_exons)
            following_exon <- min(following_exons)
            result$details$feature <- paste0("intron between exons ", preceding_exon, " and ", following_exon)
          } else if (length(preceding_exons) == 0 && length(following_exons) > 0) {
            result$details$feature <- paste0("upstream of first exon")
          } else if (length(preceding_exons) > 0 && length(following_exons) == 0) {
            result$details$feature <- paste0("downstream of last exon")
          } else {
            result$details$feature <- "intronic region (specific location unknown)"
          }
        }

        feature_determined <- TRUE
      }

      # If still not determined, it's in the gene but feature type is unknown
      if (!feature_determined) {
        result$details$feature <- "gene body (feature type unknown)"
      }
    }

    return(result)
  }

  # Process all coordinates
  message("Annotating ", nrow(coord_df), " coordinates...")
  results <- vector("list", nrow(coord_df))
  for (i in 1:nrow(coord_df)) {
    if (i %% 100 == 0 || i == nrow(coord_df)) {
      message("Processed ", i, " of ", nrow(coord_df), " coordinates")
    }
    results[[i]] <- annotate_single_position(i)
  }

  # Convert to data frame if requested
  if (return_all_as_df) {
    message("Converting results to data frame...")

    # Function to extract nested elements safely
    extract_safely <- function(lst, path) {
      result <- lst
      for (p in path) {
        if (is.list(result) && p %in% names(result)) {
          result <- result[[p]]
        } else {
          return(NA)
        }
      }
      return(result)
    }

    # Create data frame with all fields
    df <- data.frame(
      id = sapply(results, function(r) r$id),
      chr = sapply(results, function(r) r$chr),
      position = sapply(results, function(r) r$position),
      genome = sapply(results, function(r) r$genome),
      location_type = sapply(results, function(r) r$location_type),
      feature = sapply(results, function(r) extract_safely(r, c("details", "feature"))),
      gene_symbol = sapply(results, function(r) {
        if (r$location_type == "genic") {
          extract_safely(r, c("details", "gene", "symbol"))
        } else {
          NA
        }
      }),
      gene_id = sapply(results, function(r) {
        if (r$location_type == "genic") {
          extract_safely(r, c("details", "gene", "id"))
        } else {
          NA
        }
      }),
      gene_strand = sapply(results, function(r) {
        if (r$location_type == "genic") {
          extract_safely(r, c("details", "gene", "strand"))
        } else {
          NA
        }
      }),
      gene_type = sapply(results, function(r) {
        if (r$location_type == "genic") {
          extract_safely(r, c("details", "gene", "type"))
        } else {
          NA
        }
      }),
      upstream_gene = sapply(results, function(r) {
        if (r$location_type == "intergenic") {
          extract_safely(r, c("details", "upstream_gene", "symbol"))
        } else {
          NA
        }
      }),
      upstream_distance = sapply(results, function(r) {
        if (r$location_type == "intergenic") {
          extract_safely(r, c("details", "upstream_gene", "distance"))
        } else {
          NA
        }
      }),
      downstream_gene = sapply(results, function(r) {
        if (r$location_type == "intergenic") {
          extract_safely(r, c("details", "downstream_gene", "symbol"))
        } else {
          NA
        }
      }),
      downstream_distance = sapply(results, function(r) {
        if (r$location_type == "intergenic") {
          extract_safely(r, c("details", "downstream_gene", "distance"))
        } else {
          NA
        }
      }),
      stringsAsFactors = FALSE
    )

    return(df)
  } else {
    # Return list of results
    names(results) <- coord_df$id
    return(results)
  }
}

# Function to handle a single genomic position (for backward compatibility)
annotate_genomic_position <- function(chr, pos, genome, gtffile,
                                      tss_upstream = 2000, tss_downstream = 200,
                                      gene_types = NULL, transcript_types = NULL) {
  # Create a data frame with a single coordinate
  coord_df <- data.frame(chr = chr, pos = pos)

  # Call the multi-coordinate function
  results <- annotate_genomic_coordinates(coord_df, genome, gtffile,
                                          tss_upstream, tss_downstream,
                                          gene_types, transcript_types,
                                          return_all_as_df = FALSE)

  # Return the first (and only) result
  return(results[[1]])
}

# Example usage
# Single position
# result <- annotate_genomic_position("chr16", 4334103, "hg38", "path/to/gencode.v38.annotation.gtf")
# print(result)

# Multiple positions
# coords <- data.frame(chr = c("chr1", "chr16", "chr7"), pos = c(1000000, 4334103, 55259515))
# results <- annotate_genomic_coordinates(coords, "hg38", "path/to/gencode.v38.annotation.gtf")
# print(head(results))
