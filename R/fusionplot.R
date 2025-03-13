# Load required libraries
library(circlize)
library(VariantAnnotation)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(GenomeInfoDb)

# Function to read gzipped VCF and generate a Circos plot with annotations
get_fusions_from_vcf <- function(vcf_file, gtf_file = NULL, genome = "hg38",
                                 title = "", thresh = 20,
                                 allowed = paste0("chr", c(1:22)),
                                 highlight = NULL, plot = FALSE,
                                 filter = c("pass", "fully_spanned", "protein_coding"),
                                 annotate = FALSE,
                                 tss_upstream = 2000, tss_downstream = 200) {
  # Check if file exists
  if (!file.exists(vcf_file)) {
    stop(paste("Error: File not found:", vcf_file))
  }

  # Check if GTF file exists if annotation is requested
  if (annotate && is.null(gtf_file)) {
    stop("Error: GTF file required for annotation. Please provide gtf_file parameter.")
  }

  if (annotate && !file.exists(gtf_file)) {
    stop(paste("Error: GTF file not found:", gtf_file))
  }

  # Read VCF file
  message("Reading VCF file...")
  vcf <- readVcf(vcf_file, genome = genome)
  vcf <- vcf[vcf@info$SVTYPE=="BND",]
  if("pass" %in% filter){
    vcf <- vcf[vcf@fixed$FILTER=="PASS"]
  }
  if("fully_spanned" %in% filter){
    vcf <- vcf[vcf@info$NotFullySpanned=="FALSE"]
  }
  if("protein_coding" %in% filter){
    vcf <- vcf[grepl("protein_coding", sapply(vcf@info$BCSQ, function(bs){ifelse(length(bs)==0, NA, bs[1])})),]
  }
  vcf <- vcf[vcf@assays@data$DP[,1]>thresh,]

  message("Processing structural variants...")
  dat <- gsub("pbsv.BND.", "", as.character(info(vcf)$MATEID))
  parts <- strsplit(dat, "-")
  sv_data <- do.call(rbind, parts)
  sv_data <- data.frame(chr1=strsplit(sv_data[,1], ":") %>% sapply("[[", 1),
                        chr2=strsplit(sv_data[,2], ":") %>% sapply("[[", 1),
                        start1=strsplit(sv_data[,1], ":") %>% sapply("[[", 2),
                        start2=strsplit(sv_data[,2], ":") %>% sapply("[[", 2))
  sv_data <- sv_data[sv_data$chr1 %in% allowed & sv_data$chr2 %in% allowed, ]
  sv_data <- sv_data %>%
    mutate(pair_key = pmin(chr1, chr2), pair_val = pmax(chr1, chr2),
           start_key = pmin(start1, start2), start_val = pmax(start1, start2)) %>%
    distinct(pair_key, pair_val, start_key, start_val, .keep_all = TRUE) %>%
    dplyr::select(-pair_key, -pair_val, -start_key, -start_val)

  # Remove NA values (only keep inter-chromosomal events)
  sv_data <- sv_data[!is.na(sv_data$chr2) & !is.na(sv_data$start2), ]

  # Convert start positions to numeric
  sv_data$start1 <- as.numeric(sv_data$start1)
  sv_data$start2 <- as.numeric(sv_data$start2)

  # Annotate the fusion breakpoints if requested
  if (annotate) {
    message("Annotating fusion breakpoints...")
    # Prepare coordinates for annotation
    coord_df <- data.frame(
      id = 1:(nrow(sv_data) * 2),
      chr = c(sv_data$chr1, sv_data$chr2),
      pos = c(sv_data$start1, sv_data$start2)
    )

    # Annotate all breakpoints using our custom function
    annotations <- annotate_genomic_coordinates(coord_df, genome, gtf_file,
                                                tss_upstream, tss_downstream)

    # Split annotations back into breakpoint1 and breakpoint2
    n <- nrow(sv_data)
    bp1_annotations <- annotations[1:n, ]
    bp2_annotations <- annotations[(n+1):(2*n), ]

    # Add annotations to sv_data
    sv_data$bp1_gene <- bp1_annotations$gene_symbol
    sv_data$bp1_feature <- bp1_annotations$feature
    sv_data$bp1_gene_strand <- bp1_annotations$gene_strand
    sv_data$bp1_type <- bp1_annotations$location_type

    sv_data$bp2_gene <- bp2_annotations$gene_symbol
    sv_data$bp2_feature <- bp2_annotations$feature
    sv_data$bp2_gene_strand <- bp2_annotations$gene_strand
    sv_data$bp2_type <- bp2_annotations$location_type

    # For intergenic regions, add nearest genes
    sv_data$bp1_upstream_gene <- bp1_annotations$upstream_gene
    sv_data$bp1_upstream_dist <- bp1_annotations$upstream_distance
    sv_data$bp1_downstream_gene <- bp1_annotations$downstream_gene
    sv_data$bp1_downstream_dist <- bp1_annotations$downstream_distance

    sv_data$bp2_upstream_gene <- bp2_annotations$upstream_gene
    sv_data$bp2_upstream_dist <- bp2_annotations$upstream_distance
    sv_data$bp2_downstream_gene <- bp2_annotations$downstream_gene
    sv_data$bp2_downstream_dist <- bp2_annotations$downstream_distance

    # Create fusion name column (Gene1--Gene2)
    sv_data$fusion_name <- paste0(
      ifelse(sv_data$bp1_type == "genic", as.character(sv_data$bp1_gene),
             paste0(as.character(sv_data$bp1_downstream_gene), "-", as.character(sv_data$bp1_upstream_gene))),
      "--",
      ifelse(sv_data$bp2_type == "genic", as.character(sv_data$bp2_gene),
             paste0(as.character(sv_data$bp2_downstream_gene), "-", as.character(sv_data$bp2_upstream_gene)))
    )

    # Create more detailed fusion description
    sv_data$fusion_details <- paste0(
      ifelse(sv_data$bp1_type == "genic",
             paste0(sv_data$bp1_gene, " (", sv_data$bp1_feature, ")"),
             paste0("Intergenic region: nearest genes ", sv_data$bp1_downstream_gene, " and ", sv_data$bp1_upstream_gene)),
      " -- ",
      ifelse(sv_data$bp2_type == "genic",
             paste0(sv_data$bp2_gene, " (", sv_data$bp2_feature, ")"),
             paste0("Intergenic region: nearest genes ", sv_data$bp2_downstream_gene, " and ", sv_data$bp2_upstream_gene))
    )
    sv_data <- get_fusion_sequences(sv_data, n_bp = 30)

  }

  if (!plot) {
    return(sv_data)
  } else {
    message("Generating Circos plot...")

    # Prepare highlight entries if specified
    if (!is.null(highlight)) {
      highlight <- sv_data[highlight, ]
      bed3 <- highlight[, c(1, 3)]
      colnames(bed3) <- c("chr", "start")
      bed3$end <- bed3$start

      bed4 <- highlight[, c(2, 4)]
      colnames(bed4) <- c("chr", "start")
      bed4$end <- bed4$start
    }

    # Prepare all breakpoints
    bed1 <- sv_data[, c(1, 3)]
    colnames(bed1) <- c("chr", "start")
    bed1$end <- bed1$start

    bed2 <- sv_data[, c(2, 4)]
    colnames(bed2) <- c("chr", "start")
    bed2$end <- bed2$start

    # Generate colors based on annotation if available
    link_colors <- rep("grey", nrow(bed1))
    if (annotate) {
      # Color links by fusion type
      # Both breakpoints in genes
      both_genic <- sv_data$bp1_type == "genic" & sv_data$bp2_type == "genic"
      # One breakpoint in gene, one intergenic
      mixed_genic <- (sv_data$bp1_type == "genic" & sv_data$bp2_type != "genic") |
        (sv_data$bp1_type != "genic" & sv_data$bp2_type == "genic")
      # Both intergenic
      both_intergenic <- sv_data$bp1_type != "genic" & sv_data$bp2_type != "genic"

      link_colors[both_genic] <- "#E41A1C80"  # Red with transparency
      link_colors[mixed_genic] <- "#377EB880"  # Blue with transparency
      link_colors[both_intergenic] <- "#4DAF4A80"  # Green with transparency
    }

    # Initialize Circos plot
    circos.initializeWithIdeogram()

    # Add all links with appropriate colors
    circos.genomicLink(bed1, bed2, col = link_colors, border = NA)

    # Add highlighted links if specified
    if (!is.null(highlight)) {
      circos.genomicLink(bed3, bed4, col = "black", lwd = 10)

      # Add labels for highlighted fusions if annotations are available
      if (annotate) {
        # This would add text labels for the highlighted fusions
        # Implementation depends on exact visualization needs
      }
    }

    # Add title
    title(title)

    # Add legend if annotated
    if (annotate) {
      legend("bottomright",
             legend = c("Gene-Gene Fusion", "Gene-Intergenic Fusion", "Intergenic-Intergenic"),
             fill = c("#E41A1C80", "#377EB880", "#4DAF4A80"),
             title = "Fusion Types",
             border = NA,
             bty = "n")
    }

    # Clear circos
    circos.clear()
  }
  return(sv_data)
}

#' Print column names from a fusion list
#'
#' Diagnostic utility to display the structure of a fusion list,
#' showing the number of entries and available columns in each fusion set.
#'
#' @param fusion_list A list of data frames containing fusion data
#' @return NULL (prints information to console)
#' @export
#' @examples
#' # Create a sample fusion list
#' fusion_list <- list(
#'   data.frame(gene1 = c("ABI1"), gene2 = c("KMT2A")),
#'   data.frame(gene_a = c("BCR"), gene_b = c("ABL1"))
#' )
#'
#' # Print column names from the fusion list
#' print_fusion_columns(fusion_list)
print_fusion_columns <- function(fusion_list) {
  cat("Examining fusion list structure...\n")
  for (i in 1:length(fusion_list)) {
    cat(paste0("Fusion set #", i, " has ", nrow(fusion_list[[i]]), " entries with columns:\n"))
    cat(paste(colnames(fusion_list[[i]]), collapse = ", "), "\n\n")
  }
}

#' Match fusions with NeoSplice database (safe version)
#'
#' A more robust version of the fusion matcher that handles missing columns
#' and data formatting issues gracefully. This function includes additional
#' error handling and automatic column name mapping.
#'
#' @param seq_input Either a file path to sequenced fusions or a data frame containing fusion data
#' @param db_file Path to the NeoSplice database file (can be gzipped)
#' @param output_file Optional path to save results as CSV (NULL to skip saving)
#' @param breakpoint_tolerance Number of base pairs of tolerance when matching
#'   breakpoint positions (default: 100000)
#' @param verbose Whether to print progress messages (default TRUE)
#' @return Data frame of fusion matches with detailed information
#' @export
#' @examples
#' # Basic usage
#' matches <- match_neosplice_fusions_safe(
#'   seq_input = sequenced_fusions,
#'   db_file = "path/to/NeoSplice_hg38_inframe_fusion.txt.gz"
#' )
#'
#' # Process a list of fusion data frames
#' fusion_list <- list(fusions1, fusions2, fusions3)
#' all_matches <- lapply(fusion_list, function(fusion) {
#'   match_neosplice_fusions_safe(
#'     seq_input = fusion,
#'     db_file = "path/to/NeoSplice_hg38_inframe_fusion.txt.gz",
#'     breakpoint_tolerance = 200000
#'   )
#' })
match_neosplice_fusions_safe <- function(seq_input, db_file, output_file = NULL,
                                         breakpoint_tolerance = 100000, verbose = TRUE) {
  # Parse input files
  if (verbose) cat("Parsing sequenced fusions...\n")
  seq_fusions <- if (is.data.frame(seq_input)) {
    seq_input
  } else {
    tryCatch(
      parse_sequenced_fusions(seq_input),
      error = function(e) {
        message("Error parsing sequenced fusions: ", e$message)
        return(data.frame())
      }
    )
  }

  if (nrow(seq_fusions) == 0) {
    warning("No sequenced fusions found or parsing failed")
    return(data.frame())
  }

  if (verbose) cat(paste("Found", nrow(seq_fusions), "sequenced fusions\n"))

  # Print debug info about columns
  if (verbose) {
    cat("Available columns in sequenced fusions:\n")
    cat(paste(colnames(seq_fusions), collapse = ", "), "\n")
  }

  # Check if we have the minimum columns needed
  required_cols <- c("gene1", "gene2", "chr1", "chr2")
  missing_cols <- required_cols[!required_cols %in% colnames(seq_fusions)]

  if (length(missing_cols) > 0) {
    # Try to map alternative column names to required columns
    col_map <- list(
      # Common alternative names for gene1
      gene1 = c("gene_1", "gene_a", "5_gene", "upstream_gene", "first_gene"),
      # Common alternative names for gene2
      gene2 = c("gene_2", "gene_b", "3_gene", "downstream_gene", "second_gene"),
      # Common alternative names for chr1
      chr1 = c("chrom1", "chromosome1", "chr_1", "5_chr", "first_chr"),
      # Common alternative names for chr2
      chr2 = c("chrom2", "chromosome2", "chr_2", "3_chr", "second_chr")
    )

    # Try to add missing columns from alternative names
    for (req_col in missing_cols) {
      if (req_col %in% names(col_map)) {
        alt_names <- col_map[[req_col]]
        found_alt <- FALSE

        for (alt_name in alt_names) {
          if (alt_name %in% colnames(seq_fusions)) {
            seq_fusions[[req_col]] <- seq_fusions[[alt_name]]
            if (verbose) cat(paste("Mapped", alt_name, "to", req_col, "\n"))
            found_alt <- TRUE
            break
          }
        }

        if (!found_alt) {
          # If no alternative found, create an NA column
          seq_fusions[[req_col]] <- NA
          warning(paste("Could not find column", req_col, "or any alternative. Created NA column."))
        }
      }
    }
  }

  if (verbose) cat("Parsing NeoSplice fusion database...\n")
  db_fusions <- tryCatch(
    parse_neosplice_database(db_file),
    error = function(e) {
      message("Error parsing database: ", e$message)
      return(data.frame())
    }
  )

  if (nrow(db_fusions) == 0) {
    warning("No database fusions found or parsing failed")
    return(data.frame())
  }

  if (verbose) cat(paste("Found", nrow(db_fusions), "database fusions\n"))

  # Initialize results data frame
  matches <- data.frame()

  # Find matches for each sequenced fusion
  if (verbose) cat("Finding matches...\n")
  match_count <- 0

  for (i in 1:nrow(seq_fusions)) {
    seq_fusion <- seq_fusions[i, ]

    # Skip if missing essential data
    if (is.null(seq_fusion$gene1) || is.null(seq_fusion$gene2) ||
        is.null(seq_fusion$chr1) || is.null(seq_fusion$chr2) ||
        is.na(seq_fusion$gene1) || is.na(seq_fusion$gene2) ||
        is.na(seq_fusion$chr1) || is.na(seq_fusion$chr2)) {
      if (verbose) cat(paste("Skipping fusion #", i, "due to missing data\n"))
      next
    }

    # 1. Match by gene names (exact match)
    gene_matches <- tryCatch({
      db_fusions %>%
        filter(gene1 == seq_fusion$gene1 & gene2 == seq_fusion$gene2)
    }, error = function(e) {
      message("Error in gene matching: ", e$message)
      return(data.frame())
    })

    # 2. Match by chromosomes if gene matching failed or found no matches
    if (nrow(gene_matches) == 0) {
      chr_matches <- tryCatch({
        db_fusions %>%
          filter(chr1 == seq_fusion$chr1 & chr2 == seq_fusion$chr2)
      }, error = function(e) {
        message("Error in chromosome matching: ", e$message)
        return(data.frame())
      })

      if (nrow(chr_matches) > 0) {
        chr_matches$match_type <- "Chromosome match"
        matches <- bind_rows(matches, chr_matches)
        match_count <- match_count + nrow(chr_matches)
      }
    } else {
      gene_matches$match_type <- "Exact gene match"
      matches <- bind_rows(matches, gene_matches)
      match_count <- match_count + nrow(gene_matches)
    }
  }

  # Report results
  if (verbose) {
    cat(paste("Found", match_count, "potential matches\n"))

    if (match_count > 0) {
      match_summary <- matches %>%
        group_by(match_type) %>%
        summarise(count = n())

      cat("\nMatches by type:\n")
      print(match_summary)
    }
  }

  # Save results if output file is provided
  if (!is.null(output_file) && nrow(matches) > 0) {
    write.csv(matches, output_file, row.names = FALSE)
    if (verbose) cat(paste("Results saved to", output_file, "\n"))
  }

  return(matches)
}

#' Process a fusion list with the NeoSplice database
#'
#' Convenience function to process a list of fusion data frames against
#' the NeoSplice database, with error handling to continue even if some
#' elements fail.
#'
#' @param fusion_list A list of data frames containing fusion data
#' @param db_file Path to the NeoSplice database file (can be gzipped)
#' @param breakpoint_tolerance Number of base pairs of tolerance when matching
#'   breakpoint positions (default: 100000)
#' @param verbose Whether to print progress messages (default TRUE)
#' @return A list of data frames containing match results for each fusion list element
#' @export
#' @examples
#' # Create a fusion list
#' fusion_list <- list(
#'   data.frame(gene1 = "ABI1", gene2 = "KMT2A", chr1 = "chr10", chr2 = "chr11"),
#'   data.frame(gene1 = "BCR", gene2 = "ABL1", chr1 = "chr22", chr2 = "chr9")
#' )
#'
#' # Process the entire list
#' all_matches <- process_fusion_list(
#'   fusion_list,
#'   db_file = "path/to/NeoSplice_hg38_inframe_fusion.txt.gz"
#' )
process_fusion_list <- function(fusion_list, db_file,
                                breakpoint_tolerance = 1, verbose = TRUE) {
  if (!is.list(fusion_list)) {
    stop("fusion_list must be a list of data frames")
  }

  if (verbose) cat(paste("Processing", length(fusion_list), "fusion sets\n"))

  # Process each element in the list
  results <- list()

  for (i in 1:length(fusion_list)) {
    if (verbose) cat(paste("\nProcessing fusion set", i, "...\n"))

    # Skip non-data frame elements
    if (!is.data.frame(fusion_list[[i]])) {
      warning(paste("Skipping element", i, "- not a data frame"))
      results[[i]] <- data.frame()
      next
    }

    # Try to process this fusion set
    results[[i]] <- tryCatch({
      match_neosplice_fusions_safe(
        seq_input = fusion_list[[i]],
        db_file = db_file,
        breakpoint_tolerance = breakpoint_tolerance,
        verbose = verbose
      )
    }, error = function(e) {
      warning(paste("Error processing fusion set", i, ":", e$message))
      return(data.frame())
    })
  }

  # Summarize overall results
  if (verbose) {
    total_matches <- sum(sapply(results, nrow))
    cat(paste("\nTotal matches across all fusion sets:", total_matches, "\n"))

    # Count matches by type across all results
    if (total_matches > 0) {
      all_matches <- bind_rows(results)
      if (nrow(all_matches) > 0 && "match_type" %in% colnames(all_matches)) {
        match_summary <- all_matches %>%
          group_by(match_type) %>%
          summarise(count = n())

        cat("\nOverall matches by type:\n")
        print(match_summary)
      }
    }
  }

  return(results)
}#' NeoSplice Fusion Database Matcher
#'
#' This module provides specialized functions to parse and match fusion genes
#' using the NeoSplice database format. It includes functions for parsing NeoSplice
#' fusion databases, matching sequenced fusions with database entries, and
#' providing position-based matching between fusion breakpoints and gene coordinates.
#'
#' @author Your Name
#' @import dplyr
#' @import data.table
#' @name neosplice_fusion_matcher
NULL

#' Parse NeoSplice fusion database
#'
#' Specifically designed parser for the NeoSplice fusion database format.
#' This function correctly handles the column naming and format specific to NeoSplice.
#'
#' @param file_path Path to the fusion database file (can be gzipped, with .gz extension)
#' @return A data frame with standardized fusion database information
#' @export
#' @examples
#' # Parse a NeoSplice database
#' db_fusions <- parse_neosplice_database("path/to/NeoSplice_hg38_inframe_fusion.txt.gz")
#'
#' # Examine the first few rows
#' head(db_fusions)
parse_neosplice_database <- function(file_path) {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }

  # Determine if file is gzipped
  is_gzipped <- grepl("\\.gz$", file_path)

  # Read the database file
  tryCatch({
    if (is_gzipped) {
      # For gzipped files, use fread with cmd to gunzip
      library(data.table)
      db <- fread(cmd = paste("gunzip -c", shQuote(file_path)), sep = "\t",
                  header = TRUE, fill = TRUE, quote = "")
    } else {
      # For regular files
      db <- read.delim(file_path, header = TRUE, stringsAsFactors = FALSE,
                       sep = "\t", quote = "")
    }

    # Create a fusion_id column from gene1 and gene2
    if ("gene1" %in% colnames(db) && "gene2" %in% colnames(db)) {
      db$fusion_id <- paste(db$gene1, db$gene2, sep = "--")
    }

    # Map column names to standard format expected by the matcher
    col_mapping <- list(
      # NeoSplice column names to standard names
      gene1 = "gene1",
      gene2 = "gene2",
      chr1 = "gene1_chro",
      chr2 = "gene2_chro",
      gene1_txt = "gene1_txt",
      gene2_txt = "gene2_txt",
      strand1 = "gene1_strand",
      strand2 = "gene2_strand"
    )

    # Create standardized column names
    for (std_col in names(col_mapping)) {
      neosplice_col <- col_mapping[[std_col]]
      if (neosplice_col %in% colnames(db)) {
        db[[std_col]] <- db[[neosplice_col]]
      }
    }

    # Extract start and end positions from range columns
    if ("gene1_txt" %in% colnames(db)) {
      # Parse ranges in format start-end
      range_parts <- strsplit(as.character(db$gene1_txt), "-")
      db$gene1_start <- sapply(range_parts, function(x) as.numeric(x[1]))
      db$gene1_end <- sapply(range_parts, function(x) as.numeric(x[2]))
    }

    if ("gene2_txt" %in% colnames(db)) {
      # Parse ranges in format start-end
      range_parts <- strsplit(as.character(db$gene2_txt), "-")
      db$gene2_start <- sapply(range_parts, function(x) as.numeric(x[1]))
      db$gene2_end <- sapply(range_parts, function(x) as.numeric(x[2]))
    }

    # Extract information from fusion-specific columns if available
    if ("gene1_exon" %in% colnames(db)) {
      range_parts <- strsplit(as.character(db$gene1_exon), "-")
      db$gene1_exon_start <- sapply(range_parts, function(x) as.numeric(x[1]))
      db$gene1_exon_end <- sapply(range_parts, function(x) as.numeric(x[2]))
    }

    if ("gene2_exon" %in% colnames(db)) {
      range_parts <- strsplit(as.character(db$gene2_exon), "-")
      db$gene2_exon_start <- sapply(range_parts, function(x) as.numeric(x[1]))
      db$gene2_exon_end <- sapply(range_parts, function(x) as.numeric(x[2]))
    }

    return(db)
  }, error = function(e) {
    message("Error parsing NeoSplice database: ", e$message)

    # Try simpler approach for troubleshooting
    if (is_gzipped) {
      con <- gzfile(file_path, "rt")
      on.exit(close(con))
      header <- readLines(con, n = 1)
      sample_lines <- readLines(con, n = 5)
    } else {
      lines <- readLines(file_path)
      header <- lines[1]
      sample_lines <- lines[2:min(6, length(lines))]
    }

    message("Header: ", header)
    message("Sample lines: ")
    for (line in sample_lines) {
      message(line)
    }

    stop("Failed to parse NeoSplice database", call. = FALSE)
  })
}

#' Parse sequenced fusion gene data
#'
#' Parses sequenced fusion gene data from either a file path or a data frame.
#' Supports multiple input formats including 8-line format or standard tab-delimited files.
#' The function attempts to extract genomic coordinates for breakpoint analysis when available.
#' Supports both plain text and gzipped files.
#'
#' @param input Either a file path (character) or a data frame containing fusion gene data
#' @return A data frame with standardized fusion data including gene names, chromosomal locations,
#'   and genomic coordinates when available
#' @export
#' @examples
#' # From file
#' seq_fusions <- parse_sequenced_fusions("path/to/sequenced_fusions.txt")
#'
#' # From gzipped file
#' gzipped_fusions <- parse_sequenced_fusions("path/to/sequenced_fusions.txt.gz")
#'
#' # From data frame
#' df <- data.frame(
#'   gene1 = c("ABI1", "BCR"),
#'   chr1 = c("chr10", "chr22"),
#'   gene2 = c("KMT2A", "ABL1"),
#'   chr2 = c("chr11", "chr9")
#' )
#' seq_fusions <- parse_sequenced_fusions(df)
parse_sequenced_fusions <- function(input) {
  # Check if input is a file path or data frame
  if (is.character(input) && length(input) == 1) {
    # Input is a file path
    if (!file.exists(input)) {
      stop("File does not exist: ", input)
    }

    # Determine if file is gzipped
    is_gzipped <- grepl("\\.gz$", input)

    # Function to read lines from file
    read_lines_from_file <- function(file_path, is_gzipped) {
      if (is_gzipped) {
        # For gzipped files, use gzfile connection
        con <- gzfile(file_path, "rt")
        on.exit(close(con))
        lines <- readLines(con)
      } else {
        # For regular files
        lines <- readLines(file_path)
      }
      return(lines)
    }

    # Read the entire file as text
    fusion_lines <- read_lines_from_file(input, is_gzipped)

    # Determine if we're using the format from the example (8 lines per fusion)
    if (length(fusion_lines) %% 8 == 0) {
      # Process the file based on 8 lines per fusion format
      fusions <- data.frame()

      for (i in seq(1, length(fusion_lines), by = 8)) {
        if (i + 7 <= length(fusion_lines)) {
          first_line <- strsplit(fusion_lines[i], "\t")[[1]]

          if (length(first_line) >= 13) {
            fusion_entry <- data.frame(
              gene1 = first_line[1],
              chr1 = first_line[2],
              strand1 = first_line[3],
              range1 = first_line[4],
              gene2 = first_line[6],
              chr2 = first_line[7],
              strand2 = ifelse(length(first_line) >= 9, first_line[9], NA),
              range2 = ifelse(length(first_line) >= 10, first_line[10], NA),
              fusion_type = first_line[13],
              fusion_id = paste(first_line[1], first_line[6], sep = "--"),
              stringsAsFactors = FALSE
            )

            fusions <- rbind(fusions, fusion_entry)
          }
        }
      }

      return(fusions)
    } else {
      # Try to parse as tab-delimited file with headers
      message("Attempting to parse as tab-delimited file...")
      tryCatch({
        # Use read.delim for plain text or read.table with gzfile for gzipped
        if (is_gzipped) {
          con <- gzfile(input, "rt")
          on.exit(close(con))
          fusions <- read.delim(con, header = TRUE, stringsAsFactors = FALSE)
        } else {
          fusions <- read.delim(input, header = TRUE, stringsAsFactors = FALSE)
        }

        # Apply standard processing to the data frame
        return(standardize_fusion_data(fusions))
      }, error = function(e) {
        stop(paste("Failed to parse sequenced fusions file:", e$message))
      })
    }
  } else if (is.data.frame(input)) {
    # Input is a data frame
    return(standardize_fusion_data(input))
  } else {
    stop("Input must be either a file path (character) or a data frame")
  }
}

#' Standardize fusion data
#'
#' Internal function to ensure fusion data has consistent column names and formats
#'
#' @param fusions Data frame containing fusion data
#' @return Standardized fusion data frame
#' @keywords internal
standardize_fusion_data <- function(fusions) {
  # Check if we have the required columns or can create them
  required_cols <- c("gene1", "chr1", "gene2", "chr2")

  # Create or rename columns as needed
  if (!all(required_cols %in% colnames(fusions))) {
    # Try to identify columns based on common naming patterns
    col_map <- list()

    for (col in colnames(fusions)) {
      if (grepl("gene.*5|gene.*first|first.*gene|upstream.*gene|5.*gene", col, ignore.case = TRUE)) {
        col_map[["gene1"]] <- col
      } else if (grepl("gene.*3|gene.*second|second.*gene|downstream.*gene|3.*gene", col, ignore.case = TRUE)) {
        col_map[["gene2"]] <- col
      } else if (grepl("chr.*5|chr.*first|first.*chr|upstream.*chr|5.*chr", col, ignore.case = TRUE)) {
        col_map[["chr1"]] <- col
      } else if (grepl("chr.*3|chr.*second|second.*chr|downstream.*chr|3.*chr", col, ignore.case = TRUE)) {
        col_map[["chr2"]] <- col
      }
    }

    # Rename or select columns
    fusions_selected <- fusions
    for (req_col in required_cols) {
      if (req_col %in% colnames(fusions)) {
        next
      } else if (req_col %in% names(col_map)) {
        fusions_selected[[req_col]] <- fusions[[col_map[[req_col]]]]
      } else {
        warning(paste("Required column", req_col, "not found. Results may be incomplete."))
        fusions_selected[[req_col]] <- NA
      }
    }

    # Add fusion_id if not present
    if (!"fusion_id" %in% colnames(fusions_selected)) {
      fusions_selected$fusion_id <- paste(fusions_selected$gene1, fusions_selected$gene2, sep = "--")
    }

    return(fusions_selected)
  }

  # Add fusion_id if not present
  if (!"fusion_id" %in% colnames(fusions)) {
    fusions$fusion_id <- paste(fusions$gene1, fusions$gene2, sep = "--")
  }

  return(fusions)
}

#' Check if a position falls within a genomic range
#'
#' Internal function to check if a given position falls within a genomic range string.
#' Supports various range formats including "start-end", "start..end", and comma-separated formats.
#'
#' @param position Numeric position to check
#' @param range_str String representation of a genomic range
#' @return Logical indicating whether the position falls within the range
#' @keywords internal
position_in_range <- function(position, range_str) {
  if (is.na(position) || is.na(range_str)) {
    return(FALSE)
  }

  # Handle different range formats
  if (grepl("-", range_str)) {
    # Format: "start-end"
    parts <- as.numeric(strsplit(range_str, "-")[[1]])
    return(position >= min(parts) && position <= max(parts))
  } else if (grepl("\\.\\.", range_str)) {
    # Format: "start..end"
    parts <- as.numeric(strsplit(range_str, "\\.\\.")[[1]])
    return(position >= min(parts) && position <= max(parts))
  } else if (grepl(",", range_str)) {
    # Format: "start,end"
    parts <- as.numeric(strsplit(range_str, ",")[[1]])
    return(position >= min(parts) && position <= max(parts))
  } else {
    # Try to extract numbers from the string
    numbers <- as.numeric(regmatches(range_str, gregexpr("[0-9]+", range_str))[[1]])
    if (length(numbers) >= 2) {
      return(position >= min(numbers) && position <= max(numbers))
    }
  }

  # If format is not recognized, return FALSE
  return(FALSE)
}

#' Find position-based matches for fusion genes
#'
#' Identifies matches between sequenced fusion breakpoints and gene coordinates in a database.
#' Considers both direct and reciprocal matches (chr1→gene1 & chr2→gene2 or chr1→gene2 & chr2→gene1).
#'
#' @param seq_fusions Data frame of sequenced fusions with chromosome and breakpoint positions
#' @param db_fusions Data frame of database fusions with gene coordinates
#' @param position_columns Named list specifying the column names for positions in the sequenced data
#'   (default: list(chr1="chr1", pos1="start1", chr2="chr2", pos2="start2"))
#' @param gene_columns Named list specifying the column names for gene positions in the database
#'   (default: list(gene1_chr="gene1_chro", gene1_range="gene1_txt", gene2_chr="gene2_chro", gene2_range="gene2_txt"))
#' @return Data frame of matches with details on match type
#' @export
#' @examples
#' # With default column names
#' matches <- find_position_matches(seq_fusions, db_fusions)
#'
#' # With custom column specifications
#' matches <- find_position_matches(
#'   seq_fusions,
#'   db_fusions,
#'   position_columns = list(chr1="chrom1", pos1="position1", chr2="chrom2", pos2="position2"),
#'   gene_columns = list(gene1_chr="gene1_chromosome", gene1_range="gene1_coordinates",
#'                       gene2_chr="gene2_chromosome", gene2_range="gene2_coordinates")
#' )
find_position_matches <- function(seq_fusions, db_fusions,
                                  position_columns = list(chr1="chr1", pos1="start1", chr2="chr2", pos2="start2"),
                                  gene_columns = list(gene1_chr="gene1_chro", gene1_range="gene1_txt",
                                                      gene2_chr="gene2_chro", gene2_range="gene2_txt")) {

  # Validate input data frames
  if (!is.data.frame(seq_fusions) || !is.data.frame(db_fusions)) {
    stop("Both seq_fusions and db_fusions must be data frames")
  }

  # Check that required columns exist in seq_fusions
  missing_seq_cols <- position_columns[!unlist(position_columns) %in% colnames(seq_fusions)]
  if (length(missing_seq_cols) > 0) {
    stop("Missing required columns in seq_fusions: ",
         paste(unlist(missing_seq_cols), collapse=", "))
  }

  # Check that required columns exist in db_fusions
  missing_db_cols <- gene_columns[!unlist(gene_columns) %in% colnames(db_fusions)]
  if (length(missing_db_cols) > 0) {
    stop("Missing required columns in db_fusions: ",
         paste(unlist(missing_db_cols), collapse=", "))
  }

  # Initialize results data frame
  results <- data.frame()

  # Extract column names for easier reference
  chr1_col <- position_columns$chr1
  pos1_col <- position_columns$pos1
  chr2_col <- position_columns$chr2
  pos2_col <- position_columns$pos2

  gene1_chr_col <- gene_columns$gene1_chr
  gene1_range_col <- gene_columns$gene1_range
  gene2_chr_col <- gene_columns$gene2_chr
  gene2_range_col <- gene_columns$gene2_range

  # Process each sequenced fusion
  for (i in 1:nrow(seq_fusions)) {
    seq_fusion <- seq_fusions[i, ]

    # Skip if missing essential data
    if (is.na(seq_fusion[[chr1_col]]) || is.na(seq_fusion[[pos1_col]]) ||
        is.na(seq_fusion[[chr2_col]]) || is.na(seq_fusion[[pos2_col]])) {
      next
    }

    # Get position values
    chr1 <- seq_fusion[[chr1_col]]
    pos1 <- as.numeric(seq_fusion[[pos1_col]])
    chr2 <- seq_fusion[[chr2_col]]
    pos2 <- as.numeric(seq_fusion[[pos2_col]])

    # Check for direct matches: chr1→gene1 & chr2→gene2
    direct_matches <- db_fusions %>%
      filter(
        # First breakpoint matches gene1
        !!sym(gene1_chr_col) == chr1 &
          position_in_range(pos1, !!sym(gene1_range_col)) &
          # Second breakpoint matches gene2
          !!sym(gene2_chr_col) == chr2 &
          position_in_range(pos2, !!sym(gene2_range_col))
      )

    # Check for reciprocal matches: chr1→gene2 & chr2→gene1
    reciprocal_matches <- db_fusions %>%
      filter(
        # First breakpoint matches gene2
        !!sym(gene2_chr_col) == chr1 &
          position_in_range(pos1, !!sym(gene2_range_col)) &
          # Second breakpoint matches gene1
          !!sym(gene1_chr_col) == chr2 &
          position_in_range(pos2, !!sym(gene1_range_col))
      )

    # Add match type information
    if (nrow(direct_matches) > 0) {
      direct_matches$match_type <- "Direct match (chr1→gene1, chr2→gene2)"
      direct_matches$seq_fusion_id <- i  # Use row index as ID
      direct_matches$seq_chr1 <- chr1
      direct_matches$seq_pos1 <- pos1
      direct_matches$seq_chr2 <- chr2
      direct_matches$seq_pos2 <- pos2

      results <- bind_rows(results, direct_matches)
    }

    if (nrow(reciprocal_matches) > 0) {
      reciprocal_matches$match_type <- "Reciprocal match (chr1→gene2, chr2→gene1)"
      reciprocal_matches$seq_fusion_id <- i  # Use row index as ID
      reciprocal_matches$seq_chr1 <- chr1
      reciprocal_matches$seq_pos1 <- pos1
      reciprocal_matches$seq_chr2 <- chr2
      reciprocal_matches$seq_pos2 <- pos2

      results <- bind_rows(results, reciprocal_matches)
    }
  }

  return(results)
}

#' Match sequenced fusions with NeoSplice database
#'
#' This function finds matches between sequenced fusions and the NeoSplice database
#' based on gene names, chromosomes, and position ranges.
#'
#' @param seq_input Either a file path to sequenced fusions or a data frame containing fusion data
#' @param db_file Path to the NeoSplice database file (can be gzipped)
#' @param output_file Optional path to save results as CSV (NULL to skip saving)
#' @param breakpoint_tolerance Number of base pairs of tolerance when matching
#'   breakpoint positions (default: 100000)
#' @param verbose Whether to print progress messages (default TRUE)
#' @param required_columns List of column names that must be present in seq_input (default: NULL)
#' @return Data frame of fusion matches with detailed information
#' @export
#' @examples
#' # Match sequenced fusions with NeoSplice database
#' matches <- match_neosplice_fusions(
#'   seq_input = sequenced_fusions,
#'   db_file = "path/to/NeoSplice_hg38_inframe_fusion.txt.gz",
#'   breakpoint_tolerance = 200000
#' )
#'
#' # With specific required columns
#' matches <- match_neosplice_fusions(
#'   seq_input = sequenced_fusions,
#'   db_file = "path/to/NeoSplice_hg38_inframe_fusion.txt.gz",
#'   required_columns = c("gene1", "gene2", "chr1", "chr2")
#' )
match_neosplice_fusions <- function(seq_input, db_file, output_file = NULL,
                                    breakpoint_tolerance = 100000, verbose = TRUE,
                                    required_columns = NULL) {
  # Parse input files
  if (verbose) cat("Parsing sequenced fusions...\n")
  seq_fusions <- if (is.data.frame(seq_input)) {
    seq_input
  } else {
    parse_sequenced_fusions(seq_input)
  }

  if (verbose) cat(paste("Found", nrow(seq_fusions), "sequenced fusions\n"))

  # Check for required columns
  if (!is.null(required_columns)) {
    missing_cols <- required_columns[!required_columns %in% colnames(seq_fusions)]
    if (length(missing_cols) > 0) {
      stop("Missing required columns in sequenced fusions: ",
           paste(missing_cols, collapse = ", "))
    }
  }

  # Verify at minimum we have some key columns for matching
  if (!any(c("gene1", "chr1") %in% colnames(seq_fusions))) {
    warning("Input data doesn't contain gene1 or chr1 columns. Matching may not work.")
  }

  # Print debug info about columns
  if (verbose) {
    cat("Available columns in sequenced fusions:\n")
    cat(paste(colnames(seq_fusions), collapse = ", "), "\n")
  }

  if (verbose) cat("Parsing NeoSplice fusion database...\n")
  db_fusions <- parse_neosplice_database(db_file)
  if (verbose) cat(paste("Found", nrow(db_fusions), "database fusions\n"))

  # Initialize results data frame
  matches <- data.frame()

  # Find matches for each sequenced fusion
  if (verbose) cat("Finding matches...\n")

  for (i in 1:nrow(seq_fusions)) {
    seq_fusion <- seq_fusions[i, ]
    message(paste0("Processing sequenced fusion #", i))
    # # Check if required columns exist
    has_gene1 <- "bp1_gene" %in% colnames(seq_fusion) && !is.null(seq_fusion$bp1_gene)
    has_gene2 <- "bp2_gene" %in% colnames(seq_fusion) && !is.null(seq_fusion$bp2_gene)
    has_chr1 <- "chr1" %in% colnames(seq_fusion) && !is.null(seq_fusion$chr1)
    has_chr2 <- "chr2" %in% colnames(seq_fusion) && !is.null(seq_fusion$chr2)

    # # Skip if missing essential data
    # if (!has_gene1 || !has_gene2 || !has_chr1 || !has_chr2 ||
    #     is.na(seq_fusion$gene1) || is.na(seq_fusion$gene2) ||
    #     is.na(seq_fusion$chr1) || is.na(seq_fusion$chr2)) {
    #   if (verbose) cat(paste("Skipping fusion #", i, "due to missing data\n"))
    #   next
    # }

    # 1. Match by gene names (exact match)
    gene_matches <- db_fusions %>%
      filter(gene1 == seq_fusion$bp1_gene & gene2 == seq_fusion$bp2_gene)

    # 2. Match by chromosomes and breakpoint positions if available
    position_matches <- data.frame()

    if (("start1" %in% colnames(seq_fusion) || "pos1" %in% colnames(seq_fusion)) &&
        ("start2" %in% colnames(seq_fusion) || "pos2" %in% colnames(seq_fusion))) {

      # Get breakpoint positions from sequenced fusion
      pos1 <- ifelse("start1" %in% colnames(seq_fusion), as.numeric(seq_fusion$start1), as.numeric(seq_fusion$pos1))
      pos2 <- ifelse("start2" %in% colnames(seq_fusion), as.numeric(seq_fusion$start2), as.numeric(seq_fusion$pos2))

      # Skip if positions are NA
      if (!is.na(pos1) && !is.na(pos2)) {
        # Direct matches: chr1→gene1_chro & chr2→gene2_chro
        position_matches_direct <- db_fusions %>%
          filter(
            # Chromosomes match
            chr1 == seq_fusion$chr1 &
              chr2 == seq_fusion$chr2 &
              # Position 1 is within gene1's range (with tolerance)
              (pos1 >= gene1_start - breakpoint_tolerance & pos1 <= gene1_end + breakpoint_tolerance) &
              # Position 2 is within gene2's range (with tolerance)
              (pos2 >= gene2_start - breakpoint_tolerance & pos2 <= gene2_end + breakpoint_tolerance)
          )

        # Reciprocal matches: chr1→gene2_chro & chr2→gene1_chro
        position_matches_reciprocal <- db_fusions %>%
          filter(
            # Chromosomes match (in reverse)
            chr1 == seq_fusion$chr2 &
              chr2 == seq_fusion$chr1 &
              # Position 1 is within gene2's range (with tolerance)
              (pos2 >= gene1_start - breakpoint_tolerance & pos2 <= gene1_end + breakpoint_tolerance) &
              # Position 2 is within gene1's range (with tolerance)
              (pos1 >= gene2_start - breakpoint_tolerance & pos1 <= gene2_end + breakpoint_tolerance)
          )

        # Combine direct and reciprocal matches
        if (nrow(position_matches_direct) > 0) {
          position_matches_direct$match_type <- "Position match (direct)"
          position_matches <- bind_rows(position_matches, position_matches_direct)
        }

        if (nrow(position_matches_reciprocal) > 0) {
          position_matches_reciprocal$match_type <- "Position match (reciprocal)"
          position_matches <- bind_rows(position_matches, position_matches_reciprocal)
        }
      }
    }

    # 3. Match by gene name on one side and chromosome+position on the other
    hybrid_matches <- data.frame()

    if (("start1" %in% colnames(seq_fusion) || "pos1" %in% colnames(seq_fusion)) &&
        ("start2" %in% colnames(seq_fusion) || "pos2" %in% colnames(seq_fusion))) {

      # Get breakpoint positions
      pos1 <- if ("start1" %in% colnames(seq_fusion)) as.numeric(seq_fusion$start1) else as.numeric(seq_fusion$pos1)
      pos2 <- if ("start2" %in% colnames(seq_fusion)) as.numeric(seq_fusion$start2) else as.numeric(seq_fusion$pos2)

      if (!is.na(pos1) && !is.na(pos2)) {
        # Gene1 match + position2 match
        hybrid_matches_1 <- db_fusions %>%
          filter(
            gene1 == seq_fusion$bp1_gene &
              chr2 == seq_fusion$chr2 &
              (pos2 >= gene2_start - breakpoint_tolerance & pos2 <= gene2_end + breakpoint_tolerance)
          )

        # Gene2 match + position1 match
        hybrid_matches_2 <- db_fusions %>%
          filter(
            gene2 == seq_fusion$bp2_gene &
              chr1 == seq_fusion$chr1 &
              (pos1 >= gene1_start - breakpoint_tolerance & pos1 <= gene1_end + breakpoint_tolerance)
          )

        if (nrow(hybrid_matches_1) > 0) {
          hybrid_matches_1$match_type <- "Hybrid match (gene1 + position2)"
          hybrid_matches <- bind_rows(hybrid_matches, hybrid_matches_1)
        }

        if (nrow(hybrid_matches_2) > 0) {
          hybrid_matches_2$match_type <- "Hybrid match (gene2 + position1)"
          hybrid_matches <- bind_rows(hybrid_matches, hybrid_matches_2)
        }
      }
    }

    # Add match information to all matches
    all_matches <- bind_rows(
      if (nrow(gene_matches) > 0) {
        gene_matches$match_type <- "Exact gene match"
        gene_matches
      },
      position_matches,
      hybrid_matches
    )

    if (nrow(all_matches) > 0) {
      all_matches$seq_fusion_id <- i
      all_matches$seq_gene1 <- seq_fusion$gene1
      all_matches$seq_gene2 <- seq_fusion$gene2
      all_matches$seq_chr1 <- seq_fusion$chr1
      all_matches$seq_chr2 <- seq_fusion$chr2

      if ("start1" %in% colnames(seq_fusion)) all_matches$seq_pos1 <- seq_fusion$start1
      else if ("pos1" %in% colnames(seq_fusion)) all_matches$seq_pos1 <- seq_fusion$pos1

      if ("start2" %in% colnames(seq_fusion)) all_matches$seq_pos2 <- seq_fusion$start2
      else if ("pos2" %in% colnames(seq_fusion)) all_matches$seq_pos2 <- seq_fusion$pos2

      matches <- bind_rows(matches, all_matches)
    }
  }

  # Generate summary
  if (verbose) {
    cat(paste("Found", nrow(matches), "potential matches\n"))

    if (nrow(matches) > 0) {
      match_summary <- matches %>%
        group_by(match_type) %>%
        summarise(count = n())

      cat("\nMatches by type:\n")
      print(match_summary)
    }
  }

  # Save results if output file is provided
  if (!is.null(output_file) && nrow(matches) > 0) {
    write.csv(matches, output_file, row.names = FALSE)
    if (verbose) cat(paste("Results saved to", output_file, "\n"))
  }

  return(matches)
}

#' Match fusion genes by position
#'
#' Main function to identify matches between sequenced fusion breakpoints and gene coordinates
#' in a database. This function handles the entire workflow from parsing input files
#' to generating match results based on positional constraints.
#'
#' @param seq_input Either a file path to sequenced fusions or a data frame containing fusion data
#'   with chromosome and position information
#' @param db_file File path to the fusion database (can be gzipped, with .gz extension)
#' @param output_file Optional path to save results as CSV (NULL to skip saving)
#' @param position_columns Named list specifying the column names for positions in the sequenced data
#'   (default: list(chr1="chr1", pos1="start1", chr2="chr2", pos2="start2"))
#' @param gene_columns Named list specifying the column names for gene positions in the database
#'   (default: list(gene1_chr="gene1_chro", gene1_range="gene1_txt", gene2_chr="gene2_chro", gene2_range="gene2_txt"))
#' @param verbose Whether to print progress messages (default TRUE)
#' @return Data frame of fusion matches with detailed position-based match information
#' @export
#' @examples
#' # Basic usage with default column names
#' matches <- match_fusion_positions(
#'   seq_input = "sequenced_fusions.txt",
#'   db_file = "fusion_database.txt"
#' )
#'
#' # With custom column specifications
#' matches <- match_fusion_positions(
#'   seq_input = "sequenced_fusions.txt",
#'   db_file = "fusion_database.txt",
#'   position_columns = list(chr1="chrom1", pos1="position1", chr2="chrom2", pos2="position2"),
#'   gene_columns = list(gene1_chr="gene1_chromosome", gene1_range="gene1_coordinates",
#'                       gene2_chr="gene2_chromosome", gene2_range="gene2_coordinates")
#' )
match_fusion_positions <- function(seq_input, db_file, output_file = NULL,
                                   position_columns = list(chr1="chr1", pos1="start1", chr2="chr2", pos2="start2"),
                                   gene_columns = list(gene1_chr="gene1_chro", gene1_range="gene1_txt",
                                                       gene2_chr="gene2_chro", gene2_range="gene2_txt"),
                                   verbose = TRUE) {
  # Parse input files
  if (verbose) cat("Parsing sequenced breakpoint data...\n")
  seq_fusions <- if (is.data.frame(seq_input)) {
    seq_input
  } else {
    parse_sequenced_fusions(seq_input)
  }

  if (verbose) cat(paste("Found", nrow(seq_fusions), "sequenced fusion breakpoints\n"))

  if (verbose) cat("Parsing fusion database file...\n")
  db_fusions <- parse_neosplice_database(db_file)
  if (verbose) cat(paste("Found", nrow(db_fusions), "database fusion entries\n"))

  # Verify that required columns exist
  seq_cols_exist <- all(unlist(position_columns) %in% colnames(seq_fusions))
  db_cols_exist <- all(unlist(gene_columns) %in% colnames(db_fusions))

  if (!seq_cols_exist) {
    missing <- position_columns[!unlist(position_columns) %in% colnames(seq_fusions)]
    stop("Missing required columns in sequenced data: ",
         paste(unlist(missing), collapse=", "))
  }

  if (!db_cols_exist) {
    missing <- gene_columns[!unlist(gene_columns) %in% colnames(db_fusions)]
    stop("Missing required columns in database: ",
         paste(unlist(missing), collapse=", "))
  }

  # Find matches
  if (verbose) cat("Finding position-based matches...\n")
  matches <- find_position_matches(seq_fusions, db_fusions,
                                   position_columns, gene_columns)

  if (verbose) cat(paste("Found", nrow(matches), "potential position-based matches\n"))

  # Generate summary
  if (verbose && nrow(matches) > 0) {
    match_summary <- matches %>%
      group_by(match_type) %>%
      summarise(count = n())

    cat("\nSummary of matches by type:\n")
    print(match_summary)
  }

  # Save results if output file is provided
  if (!is.null(output_file) && nrow(matches) > 0) {
    write.csv(matches, output_file, row.names = FALSE)
    if (verbose) cat(paste("Results saved to", output_file, "\n"))
  }

  # Return the matches
  return(matches)
}

#' Generate Predicted Genomic Fusion Sequences
#'
#' @description
#' This function takes a data frame containing genomic fusion information and
#' extracts predicted fusion sequences from the GRCh38 reference genome. The function
#' considers the strand orientation of each gene when extracting sequences.
#'
#' @param fusion_df A data frame containing genomic fusion information with the following columns:
#'   \itemize{
#'     \item chr1: Chromosome of first breakpoint (character, e.g., "chr1")
#'     \item chr2: Chromosome of second breakpoint (character, e.g., "chr1")
#'     \item start1: Genomic position of first breakpoint (integer)
#'     \item start2: Genomic position of second breakpoint (integer)
#'     \item bp1_gene_strand: Strand of first gene ("+" or "-")
#'     \item bp2_gene_strand: Strand of second gene ("+" or "-")
#'   }
#' @param n_bp Integer specifying the number of base pairs to extract from each side
#'   of the fusion. Default is 100bp.
#'
#' @return The original data frame with an additional column 'fusion_sequence' containing
#'   the predicted fusion sequences.
#'
#' @details
#' For genes on the positive (+) strand, the function extracts n_bp downstream of the breakpoint.
#' For genes on the negative (-) strand, the function extracts n_bp upstream of the breakpoint
#' and reverse complements it.
#'
#' @examples
#' \dontrun{
#' # Create example fusion data frame
#' fusion_df <- data.frame(
#'   chr1 = c("chr1", "chr17"),
#'   chr2 = c("chr1", "chr1"),
#'   start1 = c(92294982, 70361267),
#'   start2 = c(92080390, 113497220),
#'   bp1_gene = c("GLMN", "ENSG00000267109"),
#'   bp1_gene_strand = c("-", "-"),
#'   bp2_gene = c("BTBD8", "MAGI3"),
#'   bp2_gene_strand = c("+", "+")
#' )
#'
#' # Get fusion sequences with 100bp from each side
#' results <- get_fusion_sequences(fusion_df, n_bp = 100)
#'
#' # Save results
#' write.csv(results, "fusion_sequences.csv", row.names = FALSE)
#' }
#'
#' @importFrom BSgenome.Hsapiens.UCSC.hg38 BSgenome.Hsapiens.UCSC.hg38
#' @importFrom Biostrings getSeq
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr %>%
#'
#' @export
get_fusion_sequences <- function(fusion_df, n_bp = 100) {
  # Check for required columns
  required_cols <- c("chr1", "chr2", "start1", "start2", "bp1_gene_strand", "bp2_gene_strand")
  missing_cols <- setdiff(required_cols, colnames(fusion_df))

  if (length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Reference genome
  genome <- BSgenome.Hsapiens.UCSC.hg38

  # Initialize vector to store results
  fusion_sequences <- character(nrow(fusion_df))

  for (i in 1:nrow(fusion_df)) {
    # Extract fusion info
    chr1 <- fusion_df$chr1[i]
    chr2 <- fusion_df$chr2[i]
    start1 <- fusion_df$start1[i]
    start2 <- fusion_df$start2[i]
    strand1 <- fusion_df$bp1_gene_strand[i]
    strand2 <- fusion_df$bp2_gene_strand[i]

    # Determine sequence extraction coordinates based on strand
    # Handle first breakpoint strand (default to "+" if NA)
    if (is.na(strand1)) {
      warning(paste("Row", i, "has NA in bp1_gene_strand. Defaulting to '+' strand."))
      strand1 <- "+"
    }

    if (strand1 == "+") {
      # For + strand, take downstream of breakpoint
      range1 <- GRanges(seqnames = chr1,
                        ranges = IRanges(start = start1, end = start1 + n_bp - 1),
                        strand = "+")
      seq1 <- as.character(getSeq(genome, range1))
    } else {
      # For - strand, take upstream of breakpoint, then reverse complement
      range1 <- GRanges(seqnames = chr1,
                        ranges = IRanges(start = start1 - n_bp + 1, end = start1),
                        strand = "-")
      seq1 <- as.character(getSeq(genome, range1))
    }

    # Handle second breakpoint strand (default to "+" if NA)
    if (is.na(strand2)) {
      warning(paste("Row", i, "has NA in bp2_gene_strand. Defaulting to '+' strand."))
      strand2 <- "+"
    }

    if (strand2 == "+") {
      # For + strand, take downstream of breakpoint
      range2 <- GRanges(seqnames = chr2,
                        ranges = IRanges(start = start2, end = start2 + n_bp - 1),
                        strand = "+")
      seq2 <- as.character(getSeq(genome, range2))
    } else {
      # For - strand, take upstream of breakpoint, then reverse complement
      range2 <- GRanges(seqnames = chr2,
                        ranges = IRanges(start = start2 - n_bp + 1, end = start2),
                        strand = "-")
      seq2 <- as.character(getSeq(genome, range2))
    }

    # Combine sequences to create fusion sequence
    fusion_sequences[i] <- paste0(seq1, seq2)
  }

  # Add sequences to the original dataframe
  fusion_df$predicted_fusion_sequence <- fusion_sequences

  return(fusion_df)
}

#' Create FASTA File from Fusion Sequences
#'
#' @description
#' This function creates a FASTA format file from the fusion sequences.
#'
#' @param fusion_df A data frame containing fusion information including the
#'   fusion_sequence column (output from get_fusion_sequences function).
#' @param file_path Character string specifying the output file path. Default is "fusion_sequences.fasta".
#' @param include_genes Logical indicating whether to include gene names in the FASTA headers. Default is TRUE.
#'
#' @return None, writes FASTA file to disk.
#'
#' @examples
#' \dontrun{
#' # First get fusion sequences
#' results <- get_fusion_sequences(fusion_df, n_bp = 100)
#'
#' # Then create FASTA file
#' create_fusion_fasta(results, "my_fusion_sequences.fasta")
#' }
#'
#' @export
create_fusion_fasta <- function(fusion_df, file_path = "fusion_sequences.fasta", include_genes = TRUE) {
  # Check if fusion_sequence column exists
  if (!("predicted_fusion_sequence" %in% colnames(fusion_df))) {
    stop("fusion_df must contain a 'fusion_sequence' column. Run get_fusion_sequences first.")
  }

  # Initialize vector for FASTA lines
  fasta_lines <- character()

  # Create FASTA format
  for (i in 1:nrow(fusion_df)) {
    if (include_genes && "bp1_gene" %in% colnames(fusion_df) && "bp2_gene" %in% colnames(fusion_df)) {
      header <- paste0(">Fusion_", i, "_", fusion_df$bp1_gene[i], "_", fusion_df$bp2_gene[i])
    } else {
      header <- paste0(">Fusion_", i, "_", fusion_df$chr1[i], ":", fusion_df$start1[i], "_",
                       fusion_df$chr2[i], ":", fusion_df$start2[i])
    }
    fasta_lines <- c(fasta_lines, header, fusion_df$predicted_fusion_sequence[i])
  }

  # Write to file
  writeLines(fasta_lines, file_path)
  message(paste("Created FASTA file:", file_path))
}

# Example usage (commented out)
# # First, create your dataframe (copy from your input)
# fusion_df <- data.frame(
#   chr1 = c("chr1", "chr17", "chr3", "chr8", "chr22", "chr13", "chr13", "chr9", "chr15", "chr13", "chr9"),
#   chr2 = c("chr1", "chr1", "chr11", "chr11", "chr12", "chr13", "chr13", "chr13", "chr13", "chr13", "chr13"),
#   start1 = c(92294982, 70361267, 159101394, 71487000, 41134321, 46196938, 46290388, 15978021, 92469186, 46290603, 15973060),
#   start2 = c(92080390, 113497220, 10020891, 93427007, 6686095, 32344055, 32346319, 32354081, 32595009, 42578764, 42578773),
#   bp1_gene = c("GLMN", "ENSG00000267109", "IQCJ-SCHIP1", "EYA1", "EP300", "LCP1", "LINC00563", "CCDC171", "ENSG00000309186", "LINC00563", "CCDC171"),
#   bp1_feature = c("intron", "intron", "intron", "intron", "intron", "intron", "intron", "intron", "intron", "intron", "exon 26"),
#   bp1_gene_strand = c("-", "-", "+", "-", "+", "-", "-", "+", "-", "-", "+"),
#   bp1_type = rep("genic", 11),
#   bp2_gene = c("BTBD8", "MAGI3", "SBF2", "DEUP1", "ZNF384", "BRCA2", "BRCA2", "BRCA2", "PDS5B", "TNFSF11", "TNFSF11"),
#   bp2_feature = c("exon 1", "intron", "intron", "intron", "intron", "intron", "intron", "intron", "intron", "intron", "intron"),
#   bp2_gene_strand = c("+", "+", "-", "+", "-", "+", "+", "+", "+", "+", "+")
# )
#
# # Get fusion sequences (default 100bp on each side)
# results <- get_fusion_sequences(fusion_df, n_bp = 100)
#
# # Save results to a CSV file
# write.csv(results, "fusion_sequences.csv", row.names = FALSE)
#
# # Create a FASTA file with fusion sequences
# create_fusion_fasta(results, "fusion_sequences.fasta")

# Add a function for sequence similarity calculation using string distance when Biostrings isn't available
#' Calculate sequence similarity between two sequences
#'
#' @description
#' Calculates the similarity between two nucleotide sequences using pairwise alignment
#' or string distance metrics as a fallback.
#'
#' @param seq1 First nucleotide sequence as character string
#' @param seq2 Second nucleotide sequence as character string
#' @param type Type of alignment to perform: "global" (default), "local", or "overlap"
#' @param score_only Return only the similarity score instead of detailed alignment
#'
#' @return A numeric similarity score (percent identity) between 0 and 100
#'
calculate_sequence_similarity <- function(seq1, seq2, type = "global", score_only = TRUE) {
  # Check inputs
  if (is.null(seq1) || is.null(seq2) || nchar(seq1) == 0 || nchar(seq2) == 0) {
    warning("One of the sequences is empty or NULL")
    return(NA)
  }

  # Clean sequences - remove invalid characters
  seq1 <- gsub("[^ACGTN]", "", toupper(seq1))
  seq2 <- gsub("[^ACGTN]", "", toupper(seq2))

  if (nchar(seq1) == 0 || nchar(seq2) == 0) {
    warning("After cleaning, one of the sequences is empty")
    return(NA)
  }

  # For very long sequences, use a simpler approach
  if (nchar(seq1) > 1000 || nchar(seq2) > 1000) {
    # Use string similarity metrics instead of full alignment
    # Get substrings for comparison (first 100 characters)
    sub1 <- substr(seq1, 1, min(100, nchar(seq1)))
    sub2 <- substr(seq2, 1, min(100, nchar(seq2)))

    # Calculate edit distance
    ed <- adist(sub1, sub2)[1]
    max_len <- max(nchar(sub1), nchar(sub2))
    similarity <- 100 * (1 - ed / max_len)
    return(similarity)
  }

  # Try to load Biostrings without error
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    warning("Biostrings package not available, using simpler method")
    # Use string distance as fallback
    ed <- adist(seq1, seq2)[1]
    max_len <- max(nchar(seq1), nchar(seq2))
    similarity <- 100 * (1 - ed / max_len)
    return(similarity)
  }

  # Convert to DNAString objects
  seq1_dna <- try(Biostrings::DNAString(seq1), silent = TRUE)
  seq2_dna <- try(Biostrings::DNAString(seq2), silent = TRUE)

  if (inherits(seq1_dna, "try-error") || inherits(seq2_dna, "try-error")) {
    warning("Invalid DNA sequence detected, using simpler method")
    # Use string distance as fallback
    ed <- adist(seq1, seq2)[1]
    max_len <- max(nchar(seq1), nchar(seq2))
    similarity <- 100 * (1 - ed / max_len)
    return(similarity)
  }

  # Perform pairwise alignment
  alignment <- try(
    Biostrings::pairwiseAlignment(seq1_dna, seq2_dna, type = type),
    silent = TRUE
  )

  if (inherits(alignment, "try-error")) {
    warning("Alignment failed, using simpler method")
    # Use string distance as fallback
    ed <- adist(seq1, seq2)[1]
    max_len <- max(nchar(seq1), nchar(seq2))
    similarity <- 100 * (1 - ed / max_len)
    return(similarity)
  }

  # Return score or full alignment
  if (score_only) {
    return(Biostrings::pid(alignment, type = "PID1"))  # Percent identity
  } else {
    return(alignment)
  }
}

#' Match sequenced fusions with NeoSplice database using sequence similarity
#'
#' Enhanced version of match_neosplice_fusions that adds sequence similarity comparison
#' between predicted fusion sequences and database fusion sequences.
#'
#' @param seq_input Either a file path to sequenced fusions or a data frame containing fusion data
#' @param db_file Path to the NeoSplice database file (can be gzipped)
#' @param output_file Optional path to save results as CSV (NULL to skip saving)
#' @param breakpoint_tolerance Number of base pairs of tolerance when matching
#'   breakpoint positions (default: 100000)
#' @param sequence_similarity_threshold Minimum percent identity required for sequence match (default: 60)
#' @param verbose Whether to print progress messages (default TRUE)
#' @param required_columns List of column names that must be present in seq_input (default: NULL)
#' @param add_sequence_match Whether to find the best matching sequences (default: TRUE)
#' @param seq_column Name of the column containing predicted fusion sequences (default: "predicted_fusion_sequence")
#' @param max_candidates Maximum number of candidates to test for sequence similarity (default: 50)
#'
#' @return Data frame of fusion matches with detailed information including sequence matches
#' @export
#'
match_neosplice_fusions_with_sequence <- function(seq_input, db_file, output_file = NULL,
                                                  breakpoint_tolerance = 100000,
                                                  sequence_similarity_threshold = 60,
                                                  verbose = TRUE,
                                                  required_columns = NULL,
                                                  add_sequence_match = TRUE,
                                                  seq_column = "predicted_fusion_sequence",
                                                  max_candidates = 50) {
  # Parse input files
  if (verbose) cat("Parsing sequenced fusions...\n")
  seq_fusions <- if (is.data.frame(seq_input)) {
    seq_input
  } else {
    parse_sequenced_fusions(seq_input)
  }

  if (verbose) cat(paste("Found", nrow(seq_fusions), "sequenced fusions\n"))

  # Check for required columns
  if (!is.null(required_columns)) {
    missing_cols <- required_columns[!required_columns %in% colnames(seq_fusions)]
    if (length(missing_cols) > 0) {
      stop("Missing required columns in sequenced fusions: ",
           paste(missing_cols, collapse = ", "))
    }
  }

  # Print debug info about columns
  if (verbose) {
    cat("Available columns in sequenced fusions:\n")
    cat(paste(colnames(seq_fusions), collapse = ", "), "\n")
  }

  # Check if sequence column exists
  if (add_sequence_match && !seq_column %in% colnames(seq_fusions)) {
    warning(paste("Sequence column", seq_column, "not found in input data. Sequence matching will be skipped."))
    add_sequence_match <- FALSE
  }

  if (verbose) cat("Parsing NeoSplice fusion database...\n")
  db_fusions <- parse_neosplice_database(db_file)
  if (verbose) cat(paste("Found", nrow(db_fusions), "database fusions\n"))

  # Check if fuse_contig column exists in database
  if (add_sequence_match && !"fuse_contig" %in% colnames(db_fusions)) {
    warning("'fuse_contig' column not found in database. Available columns:")
    cat(paste(colnames(db_fusions), collapse = ", "))
    add_sequence_match <- FALSE
  }

  # Initialize results data frame
  matches <- data.frame()

  # Find matches for each sequenced fusion
  if (verbose) cat("Finding matches...\n")

  for (i in 1:nrow(seq_fusions)) {
    seq_fusion <- seq_fusions[i, ]
    if (verbose) cat(paste0("Processing sequenced fusion #", i, "\n"))

    # Check if required columns exist
    has_gene1 <- "bp1_gene" %in% colnames(seq_fusion) && !is.na(seq_fusion$bp1_gene)
    has_gene2 <- "bp2_gene" %in% colnames(seq_fusion) && !is.na(seq_fusion$bp2_gene)
    has_chr1 <- "chr1" %in% colnames(seq_fusion) && !is.na(seq_fusion$chr1)
    has_chr2 <- "chr2" %in% colnames(seq_fusion) && !is.na(seq_fusion$chr2)

    # 1. Match by gene names (exact match)
    gene_matches <- data.frame()
    if (has_gene1 && has_gene2) {
      gene_matches <- tryCatch({
        db_fusions %>%
          filter(gene1 == seq_fusion$bp1_gene & gene2 == seq_fusion$bp2_gene)
      }, error = function(e) {
        message("Error in gene matching: ", e$message)
        return(data.frame())
      })

      if (nrow(gene_matches) > 0) {
        gene_matches$match_type <- "Exact gene match"
      }
    }

    # 2. Match by chromosomes and breakpoint positions
    pos_matches <- data.frame()
    if (has_chr1 && has_chr2 && "start1" %in% colnames(seq_fusion) && "start2" %in% colnames(seq_fusion)) {
      pos1 <- as.numeric(seq_fusion$start1)
      pos2 <- as.numeric(seq_fusion$start2)

      # Only proceed if positions are not NA
      if (!is.na(pos1) && !is.na(pos2)) {
        pos_matches <- tryCatch({
          # Direct matches: chr1→gene1_chro & chr2→gene2_chro
          db_fusions %>%
            filter(
              # Chromosomes match
              gene1_chro == seq_fusion$chr1 &
                gene2_chro == seq_fusion$chr2 &
                # Positions are near gene locations
                (pos1 >= as.numeric(sub("-.*$", "", gene1_txt)) - breakpoint_tolerance &
                   pos1 <= as.numeric(sub("^.*-", "", gene1_txt)) + breakpoint_tolerance) &
                (pos2 >= as.numeric(sub("-.*$", "", gene2_txt)) - breakpoint_tolerance &
                   pos2 <= as.numeric(sub("^.*-", "", gene2_txt)) + breakpoint_tolerance)
            )
        }, error = function(e) {
          message("Error in position matching: ", e$message)
          return(data.frame())
        })

        if (nrow(pos_matches) > 0) {
          pos_matches$match_type <- "Position match"
        }
      }
    }

    # 3. Sequence matching if enabled
    seq_matches <- data.frame()
    if (add_sequence_match && seq_column %in% colnames(seq_fusion)) {
      # Get the predicted fusion sequence
      pred_seq <- seq_fusion[[seq_column]]

      if (verbose) cat(paste("  Row", i, "sequence available:", !is.na(pred_seq), ", length:", ifelse(is.na(pred_seq), 0, nchar(pred_seq)), "\n"))

      if (!is.na(pred_seq) && nchar(pred_seq) > 0) {
        # First try to match with gene matches for efficiency
        candidates <- if (nrow(gene_matches) > 0) {
          if (verbose) cat("  Using gene matches as candidates\n")
          gene_matches
        } else if (nrow(pos_matches) > 0) {
          if (verbose) cat("  Using position matches as candidates\n")
          pos_matches
        } else {
          if (verbose) cat("  No matches yet, filtering database...\n")
          # If no other matches, try the entire database (limited for performance)
          # Use gene/chr filter to narrow down candidates
          filtered_db <- db_fusions
          if (has_gene1) {
            filtered_db <- filtered_db %>% filter(gene1 == seq_fusion$bp1_gene | gene2 == seq_fusion$bp1_gene)
          } else if (has_gene2) {
            filtered_db <- filtered_db %>% filter(gene1 == seq_fusion$bp2_gene | gene2 == seq_fusion$bp2_gene)
          } else if (has_chr1) {
            filtered_db <- filtered_db %>% filter(gene1_chro == seq_fusion$chr1 | gene2_chro == seq_fusion$chr1)
          } else if (has_chr2) {
            filtered_db <- filtered_db %>% filter(gene1_chro == seq_fusion$chr2 | gene2_chro == seq_fusion$chr2)
          }

          # Take a subset if too many candidates
          if (nrow(filtered_db) > max_candidates) {
            if (verbose) cat(paste("  Sampling from", nrow(filtered_db), "candidates\n"))
            filtered_db <- filtered_db[sample(nrow(filtered_db), max_candidates), ]
          } else {
            if (verbose) cat(paste("  Using all", nrow(filtered_db), "filtered candidates\n"))
          }
          filtered_db
        }

        # Check if fuse_contig column exists
        if (!"fuse_contig" %in% colnames(candidates)) {
          if (verbose) cat("  'fuse_contig' column not found in database. Available columns:\n")
          if (verbose) cat(paste("  ", paste(colnames(candidates), collapse=", "), "\n"))
        } else {
          if (verbose) cat(paste("  Found", nrow(candidates), "candidates with fuse_contig column\n"))

          # Calculate sequence similarity for candidates
          if (nrow(candidates) > 0) {
            # Calculate similarity scores for each candidate
            similarity_scores <- numeric(nrow(candidates))
            success_count <- 0
            error_count <- 0

            for (j in 1:nrow(candidates)) {
              db_seq <- candidates$fuse_contig[j]
              if (!is.na(db_seq) && nchar(db_seq) > 0) {
                # Use tryCatch to avoid stopping on errors
                similarity_scores[j] <- tryCatch({
                  score <- calculate_sequence_similarity(pred_seq, db_seq)
                  success_count <- success_count + 1
                  score
                }, error = function(e) {
                  error_count <- error_count + 1
                  if (verbose && error_count <= 3) cat(paste("  Error calculating similarity:", e$message, "\n"))
                  NA
                })
              } else {
                similarity_scores[j] <- NA
              }
            }

            if (verbose) cat(paste("  Calculated", success_count, "similarity scores,", error_count, "errors\n"))

            # Add similarity scores to candidates
            candidates$sequence_similarity <- similarity_scores

            # Filter candidates by similarity threshold
            seq_matches <- candidates %>%
              filter(!is.na(sequence_similarity) & sequence_similarity >= sequence_similarity_threshold) %>%
              arrange(desc(sequence_similarity))

            if (nrow(seq_matches) > 0) {
              if (verbose) cat(paste("  Found", nrow(seq_matches), "sequence matches above threshold\n"))
              # Add match type
              if ("match_type" %in% colnames(seq_matches)) {
                seq_matches$match_type <- paste0(seq_matches$match_type, " + Sequence match")
              } else {
                seq_matches$match_type <- "Sequence match"
              }
            } else {
              if (verbose) cat("  No sequence matches above threshold\n")
            }
          }
        }
      }
    }

    # Combine all match types
    all_matches <- bind_rows(
      if (nrow(gene_matches) > 0) gene_matches,
      if (nrow(pos_matches) > 0) pos_matches,
      if (nrow(seq_matches) > 0) seq_matches
    ) %>% distinct()

    if (nrow(all_matches) > 0) {
      # Add sequenced fusion info to matches
      all_matches$seq_fusion_id <- i
      all_matches$seq_gene1 <- seq_fusion$bp1_gene
      all_matches$seq_gene2 <- seq_fusion$bp2_gene
      all_matches$seq_chr1 <- seq_fusion$chr1
      all_matches$seq_chr2 <- seq_fusion$chr2
      all_matches$seq_pos1 <- seq_fusion$start1
      all_matches$seq_pos2 <- seq_fusion$start2

      # Add to overall matches
      matches <- bind_rows(matches, all_matches)
    }
  }

  # Generate summary
  if (verbose) {
    cat(paste("Found", nrow(matches), "potential matches\n"))

    if (nrow(matches) > 0) {
      match_summary <- matches %>%
        group_by(match_type) %>%
        summarise(count = n())

      cat("\nMatches by type:\n")
      print(match_summary)
    }
  }

  # Save results if output file is provided
  if (!is.null(output_file) && nrow(matches) > 0) {
    write.csv(matches, output_file, row.names = FALSE)
    if (verbose) cat(paste("Results saved to", output_file, "\n"))
  }

  return(matches)
}

# Usage example:
# match <- match_neosplice_fusions_with_sequence(
#   seq_input = fusions_S1,
#   db_file = "path/to/NeoSplice_hg38_inframe_fusion.txt.gz",
#   breakpoint_tolerance = 2,
#   sequence_similarity_threshold = 60,
#   add_sequence_match = TRUE,
#   seq_column = "predicted_fusion_sequence"
# )
