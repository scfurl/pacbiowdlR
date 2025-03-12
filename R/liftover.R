#' Lift over fusion gene coordinates from hg19 to hg38
#'
#' This function takes a fusion gene dataset in the ChimerDB4 format with hg19 coordinates
#' and lifts over the coordinates to hg38 using the UCSC chain file.
#'
#' @param input_file Character string specifying the path to the input file. File should be
#'   tab-delimited in ChimerDB4 format.
#' @param output_dir Character string specifying the directory where output files will be written.
#'   Defaults to the current working directory.
#' @param chain_file Character string specifying the path to the chain file. If NULL (default),
#'   the function will download the chain file from UCSC.
#' @param column_names Character vector of column names for the input file. If NULL (default),
#'   the function will use the default ChimerDB4 column names.
#'
#' @return A list containing three data frames:
#' \describe{
#'   \item{combined}{Data frame with both hg19 and hg38 coordinates}
#'   \item{hg38_only}{Data frame with only hg38 coordinates for successful liftovers}
#'   \item{failed}{Data frame with records that failed to lift over}
#' }
#'
#' @details The function creates three output files in the specified directory:
#' \describe{
#'   \item{chimer_data_hg19_and_hg38.txt}{Contains both hg19 and hg38 coordinates}
#'   \item{chimer_data_hg38.txt}{Contains only hg38 coordinates for successful liftovers}
#'   \item{failed_liftover.txt}{Contains records that failed to lift over}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic usage with default parameters
#' liftover_fusion_coordinates("path/to/chimerdb4_file.txt")
#'
#' # Specify custom output directory and column names
#' result <- liftover_fusion_coordinates(
#'   "path/to/chimerdb4_file.txt",
#'   output_dir = "path/to/output",
#'   column_names = c("DB", "Cancer", "Sample", "Gene1", "Chr1", "Pos1", "Strand1",
#'                    "Gene2", "Chr2", "Pos2", "Strand2")
#' )
#'
#' # Use a local chain file
#' liftover_fusion_coordinates(
#'   "path/to/chimerdb4_file.txt",
#'   chain_file = "path/to/hg19ToHg38.over.chain"
#' )
#' }
#'
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr mutate left_join filter select n
#' @importFrom readr read_delim write_delim
#' @importFrom utils download.file
#'
#' @export
liftover_fusion_coordinates <- function(input_file,
                                        output_dir = ".",
                                        chain_file = NULL,
                                        column_names = NULL) {
  # Check if input file exists
  if (!file.exists(input_file)) {
    stop("Input file does not exist: ", input_file)
  }

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Default column names for ChimerDB4 format
  if (is.null(column_names)) {
    column_names <- c("DB", "Cancer_Type", "Sample", "Gene1", "Chr1", "Pos1", "Strand1",
                      "Gene2", "Chr2", "Pos2", "Strand2")
  }

  # Read in the ChimerDB4 data
  chimer_data <- read.table(input_file, sep="\t", header=FALSE, stringsAsFactors=FALSE)

  # Assign column names
  # If there are more columns in the file than in column_names, extend column_names
  if (ncol(chimer_data) > length(column_names)) {
    extra_cols <- ncol(chimer_data) - length(column_names)
    column_names <- c(column_names, paste0("V", (length(column_names) + 1):(length(column_names) + extra_cols)))
  }
  # If there are fewer columns in the file than in column_names, truncate column_names
  if (ncol(chimer_data) < length(column_names)) {
    column_names <- column_names[1:ncol(chimer_data)]
  }

  colnames(chimer_data) <- column_names

  # Ensure the essential columns exist
  required_cols <- c("Chr1", "Pos1", "Strand1", "Chr2", "Pos2", "Strand2")
  missing_cols <- required_cols[!required_cols %in% colnames(chimer_data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Normalize strand values to ensure they are valid ('+', '-', or '*')
  normalize_strand <- function(strand_values) {
    # Convert to character to handle factors safely
    strand_values <- as.character(strand_values)

    # Replace invalid values with '*' (any strand)
    valid_values <- c("+", "-", "*")
    strand_values[!strand_values %in% valid_values] <- "*"

    return(strand_values)
  }

  # Normalize strand values
  chimer_data$Strand1 <- normalize_strand(chimer_data$Strand1)
  chimer_data$Strand2 <- normalize_strand(chimer_data$Strand2)

  # Create two GRanges objects - one for each fusion partner
  gr1 <- GRanges(
    seqnames = chimer_data$Chr1,
    ranges = IRanges(start = as.numeric(chimer_data$Pos1), width = 1),
    strand = chimer_data$Strand1,
    fusion_id = 1:nrow(chimer_data),
    gene = if ("Gene1" %in% colnames(chimer_data)) chimer_data$Gene1 else NA
  )

  gr2 <- GRanges(
    seqnames = chimer_data$Chr2,
    ranges = IRanges(start = as.numeric(chimer_data$Pos2), width = 1),
    strand = chimer_data$Strand2,
    fusion_id = 1:nrow(chimer_data),
    gene = if ("Gene2" %in% colnames(chimer_data)) chimer_data$Gene2 else NA
  )

  # Handle the chain file
  if (is.null(chain_file)) {
    chain_file <- file.path(output_dir, "hg19ToHg38.over.chain")
    if (!file.exists(chain_file)) {
      temp_gz <- file.path(output_dir, "hg19ToHg38.over.chain.gz")
      download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                    temp_gz)
      system(paste("gunzip", temp_gz))
    }
  } else if (!file.exists(chain_file)) {
    stop("Chain file does not exist: ", chain_file)
  }

  # Import the chain file
  chain <- import.chain(chain_file)

  # Perform liftOver
  gr1_hg38 <- liftOver(gr1, chain)
  gr2_hg38 <- liftOver(gr2, chain)

  # Convert liftOver results to GRanges objects
  # Some coordinates might not lift over, so we need to handle that
  gr1_hg38_unlisted <- unlist(gr1_hg38)
  gr2_hg38_unlisted <- unlist(gr2_hg38)

  # Create dataframes with the lifted coordinates
  if (length(gr1_hg38_unlisted) > 0) {
    df1_hg38 <- data.frame(
      fusion_id = gr1_hg38_unlisted$fusion_id,
      chr1_hg38 = seqnames(gr1_hg38_unlisted),
      pos1_hg38 = start(gr1_hg38_unlisted),
      strand1_hg38 = strand(gr1_hg38_unlisted),
      gene1 = gr1_hg38_unlisted$gene
    )
  } else {
    df1_hg38 <- data.frame(
      fusion_id = integer(),
      chr1_hg38 = character(),
      pos1_hg38 = integer(),
      strand1_hg38 = character(),
      gene1 = character()
    )
  }

  if (length(gr2_hg38_unlisted) > 0) {
    df2_hg38 <- data.frame(
      fusion_id = gr2_hg38_unlisted$fusion_id,
      chr2_hg38 = seqnames(gr2_hg38_unlisted),
      pos2_hg38 = start(gr2_hg38_unlisted),
      strand2_hg38 = strand(gr2_hg38_unlisted),
      gene2 = gr2_hg38_unlisted$gene
    )
  } else {
    df2_hg38 <- data.frame(
      fusion_id = integer(),
      chr2_hg38 = character(),
      pos2_hg38 = integer(),
      strand2_hg38 = character(),
      gene2 = character()
    )
  }

  # Join the results with the original data
  result1 <- chimer_data %>%
    mutate(fusion_id = 1:n()) %>%
    left_join(df1_hg38, by = "fusion_id")

  result2 <- result1 %>%
    left_join(df2_hg38, by = "fusion_id")

  # Identify records where liftOver failed
  failed_liftover <- result2 %>%
    filter(is.na(chr1_hg38) | is.na(chr2_hg38))

  if (nrow(failed_liftover) > 0) {
    warning(paste0(nrow(failed_liftover), " records failed to lift over."))
    # Write failed records to a file
    write.table(failed_liftover, file.path(output_dir, "failed_liftover.txt"),
                sep="\t", quote=FALSE, row.names=FALSE)
  }

  # Create a final table with both hg19 and hg38 coordinates
  gene_cols <- intersect(c("Gene1", "Gene2"), colnames(chimer_data))
  base_cols <- setdiff(colnames(chimer_data), c("Chr1", "Pos1", "Strand1", "Chr2", "Pos2", "Strand2", gene_cols))

  final_result <- result2 %>%
    select(all_of(c(
      base_cols,
      "Gene1", "Chr1", "Pos1", "Strand1", "chr1_hg38", "pos1_hg38", "strand1_hg38",
      "Gene2", "Chr2", "Pos2", "Strand2", "chr2_hg38", "pos2_hg38", "strand2_hg38"
    )))

  # Write the result to a file
  write.table(final_result, file.path(output_dir, "chimer_data_hg19_and_hg38.txt"),
              sep="\t", quote=FALSE, row.names=FALSE)

  # Create a hg38-only version
  select_cols <- c(base_cols)

  if ("Gene1" %in% colnames(chimer_data)) {
    select_cols <- c(select_cols, c("gene1", "chr1_hg38", "pos1_hg38", "strand1_hg38"))
  } else {
    select_cols <- c(select_cols, c("chr1_hg38", "pos1_hg38", "strand1_hg38"))
  }

  if ("Gene2" %in% colnames(chimer_data)) {
    select_cols <- c(select_cols, c("gene2", "chr2_hg38", "pos2_hg38", "strand2_hg38"))
  } else {
    select_cols <- c(select_cols, c("chr2_hg38", "pos2_hg38", "strand2_hg38"))
  }

  hg38_only <- result2 %>%
    filter(!is.na(chr1_hg38) & !is.na(chr2_hg38)) %>%
    select(all_of(select_cols))

  # Rename columns
  names(hg38_only) <- gsub("_hg38", "", names(hg38_only))
  if ("gene1" %in% names(hg38_only)) names(hg38_only)[names(hg38_only) == "gene1"] <- "Gene1"
  if ("gene2" %in% names(hg38_only)) names(hg38_only)[names(hg38_only) == "gene2"] <- "Gene2"

  # Write the hg38-only result to a file
  write.table(hg38_only, file.path(output_dir, "chimer_data_hg38.txt"),
              sep="\t", quote=FALSE, row.names=FALSE)

  # Print summary statistics
  message("Total records: ", nrow(chimer_data))
  message("Successfully lifted over: ", nrow(hg38_only))
  message("Failed to lift over: ", nrow(chimer_data) - nrow(hg38_only))

  # Return all results as a list
  return(list(
    combined = final_result,
    hg38_only = hg38_only,
    failed = failed_liftover
  ))
}

#' Compare sequenced fusions with liftover database
#'
#' This function compares fusion breakpoints from sequencing data with fusions from a
#' liftover database (e.g., ChimerDB4) to identify potential matches. It can accept
#' both file paths or data frames as input.
#'
#' @param sequenced_fusions File path or data frame containing sequenced fusions data
#' @param liftover_db File path or data frame containing lifted over ChimerDB4 data
#' @param window_size Number of base pairs around each breakpoint to consider as a match (default: 10000)
#' @param output_file Path to write matching results (optional)
#'
#' @return A data frame containing matching fusions
#'
#' @importFrom dplyr filter mutate select left_join bind_rows
#' @importFrom readr read_delim write_csv
#' @export
compare_fusions <- function(sequenced_fusions, liftover_db, window_size = 10000, output_file = NULL) {
  # Read sequenced fusions if a file path is provided
  if (is.character(sequenced_fusions) && length(sequenced_fusions) == 1) {
    message("Reading sequenced fusions file...")
    seq_fusions <- read.table(sequenced_fusions, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else if (is.data.frame(sequenced_fusions)) {
    message("Using provided sequenced fusions data frame...")
    seq_fusions <- sequenced_fusions
  } else {
    stop("sequenced_fusions must be either a file path or a data frame")
  }

  # Read liftover database if a file path is provided
  if (is.character(liftover_db) && length(liftover_db) == 1) {
    message("Reading liftover database file...")
    liftover_data <- read.table(liftover_db, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  } else if (is.data.frame(liftover_db)) {
    message("Using provided liftover database data frame...")
    liftover_data <- liftover_db
  } else {
    stop("liftover_db must be either a file path or a data frame")
  }

  # Normalize chromosome names in both datasets
  message("Normalizing chromosome names...")
  normalize_chr <- function(chr) {
    chr <- as.character(chr)
    # Add "chr" prefix if not present
    chr <- ifelse(!grepl("^chr", chr, ignore.case = TRUE), paste0("chr", chr), chr)
    # Convert chrM to chrMT if needed
    chr <- ifelse(chr == "chrM", "chrMT", chr)
    # Standardize case
    chr <- toupper(chr)
    return(chr)
  }

  # Normalize sequenced fusion chromosome names
  if("chr1" %in% colnames(seq_fusions)) {
    seq_fusions$chr1 <- normalize_chr(seq_fusions$chr1)
    seq_fusions$chr2 <- normalize_chr(seq_fusions$chr2)
  }

  # Normalize liftover database chromosome names depending on column names
  if("Chr1" %in% colnames(liftover_data)) {
    liftover_data$Chr1 <- normalize_chr(liftover_data$Chr1)
    liftover_data$Chr2 <- normalize_chr(liftover_data$Chr2)
  } else if("chr1" %in% colnames(liftover_data)) {
    liftover_data$chr1 <- normalize_chr(liftover_data$chr1)
    liftover_data$chr2 <- normalize_chr(liftover_data$chr2)
  }

  # Determine column names for positions in liftover data
  pos1_col <- if("Pos1" %in% colnames(liftover_data)) "Pos1" else "pos1"
  pos2_col <- if("Pos2" %in% colnames(liftover_data)) "Pos2" else "pos2"

  # Determine column names for chromosomes in liftover data
  chr1_col <- if("Chr1" %in% colnames(liftover_data)) "Chr1" else "chr1"
  chr2_col <- if("Chr2" %in% colnames(liftover_data)) "Chr2" else "chr2"

  # Find matches
  message("Finding fusion matches...")
  matches <- data.frame()

  for(i in 1:nrow(seq_fusions)) {
    seq_fusion <- seq_fusions[i,]

    # Get coordinates
    seq_chr1 <- seq_fusion$chr1
    seq_chr2 <- seq_fusion$chr2
    seq_pos1 <- seq_fusion$start1
    seq_pos2 <- seq_fusion$start2

    # Find matches in database
    # First, filter by chromosomes
    chr_matches <- liftover_data[
      (liftover_data[[chr1_col]] == seq_chr1 & liftover_data[[chr2_col]] == seq_chr2) |
        (liftover_data[[chr1_col]] == seq_chr2 & liftover_data[[chr2_col]] == seq_chr1),
    ]

    if(nrow(chr_matches) > 0) {
      # Check positions within window
      position_matches <- chr_matches[
        # Direct match: chr1-chr2
        (abs(as.numeric(chr_matches[[pos1_col]]) - seq_pos1) <= window_size &
           abs(as.numeric(chr_matches[[pos2_col]]) - seq_pos2) <= window_size) |
          # Reverse match: chr2-chr1
          (abs(as.numeric(chr_matches[[pos1_col]]) - seq_pos2) <= window_size &
             abs(as.numeric(chr_matches[[pos2_col]]) - seq_pos1) <= window_size),
      ]

      if(nrow(position_matches) > 0) {
        # Add sequenced fusion information and calculate distances
        for(j in 1:nrow(position_matches)) {
          # Add sequenced fusion details
          position_matches$seq_fusion_name[j] <- if("fusion_name" %in% colnames(seq_fusion)) seq_fusion$fusion_name else paste(seq_chr1, seq_pos1, seq_chr2, seq_pos2, sep="-")
          position_matches$seq_chr1[j] <- seq_chr1
          position_matches$seq_pos1[j] <- seq_pos1
          position_matches$seq_chr2[j] <- seq_chr2
          position_matches$seq_pos2[j] <- seq_pos2

          # Add gene information if available
          if("bp1_gene" %in% colnames(seq_fusion)) position_matches$seq_gene1[j] <- seq_fusion$bp1_gene
          if("bp2_gene" %in% colnames(seq_fusion)) position_matches$seq_gene2[j] <- seq_fusion$bp2_gene

          # Calculate distances
          if(position_matches[[chr1_col]][j] == seq_chr1 && position_matches[[chr2_col]][j] == seq_chr2) {
            # Direct match
            position_matches$dist1[j] <- abs(as.numeric(position_matches[[pos1_col]][j]) - seq_pos1)
            position_matches$dist2[j] <- abs(as.numeric(position_matches[[pos2_col]][j]) - seq_pos2)
            position_matches$orientation[j] <- "direct"
          } else {
            # Reverse match
            position_matches$dist1[j] <- abs(as.numeric(position_matches[[pos1_col]][j]) - seq_pos2)
            position_matches$dist2[j] <- abs(as.numeric(position_matches[[pos2_col]][j]) - seq_pos1)
            position_matches$orientation[j] <- "reverse"
          }

          position_matches$max_dist[j] <- max(position_matches$dist1[j], position_matches$dist2[j])
          position_matches$sum_dist[j] <- position_matches$dist1[j] + position_matches$dist2[j]
        }

        # Add to matches
        matches <- rbind(matches, position_matches)
      }
    }
  }

  # Sort by distance
  if(nrow(matches) > 0) {
    matches <- matches[order(matches$max_dist),]

    # Write output if requested
    if(!is.null(output_file)) {
      write.csv(matches, output_file, row.names = FALSE)
      message("Matches written to ", output_file)
    }

    message("Found ", nrow(matches), " matching fusions within ", window_size, "bp window")
  } else {
    message("No matching fusions found within ", window_size, "bp window")
  }

  return(matches)
}

#' Run the fusion comparison with options to adjust matching parameters
#'
#' @param seq_input File path or data frame containing sequenced fusions
#' @param liftover_input File path or data frame containing liftover database
#' @param output_file Path to write results (optional)
#' @param window_sizes Vector of window sizes to test (default: c(1000, 5000, 10000, 50000))
#'
#' @return A list of match data frames for each window size
#'
#' @export
run_fusion_comparison <- function(seq_input, liftover_input, output_file = NULL,
                                  window_sizes = c(1000, 5000, 10000, 50000)) {
  results <- list()

  for(window in window_sizes) {
    message("\nChecking with window size: ", window, "bp")
    if(!is.null(output_file)) {
      window_output <- paste0(
        tools::file_path_sans_ext(output_file),
        "_window", window, ".",
        tools::file_ext(output_file)
      )
    } else {
      window_output <- NULL
    }

    results[[paste0("window_", window)]] <- compare_fusions(
      sequenced_fusions = seq_input,
      liftover_db = liftover_input,
      window_size = window,
      output_file = window_output
    )
  }

  return(results)
}

# Example usage
#
# # With file paths:
# results <- run_fusion_comparison(
#   seq_input = "path/to/sequenced_fusions.txt",
#   liftover_input = "path/to/chimer_data_hg38.txt",
#   output_file = "fusion_matches.csv"
# )
#
# # With data frames:
# seq_data <- read.table("path/to/sequenced_fusions.txt", header=TRUE, sep="\t")
# liftover_data <- read.table("path/to/chimer_data_hg38.txt", header=TRUE, sep="\t")
#
# results <- run_fusion_comparison(
#   seq_input = seq_data,
#   liftover_input = liftover_data
# )
