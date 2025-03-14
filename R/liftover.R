#' Lift over fusion gene coordinates from hg19 to hg38 (multi-file version)
#'
#' This function accepts *multiple* input files (a character vector) in ChimerDB4-like format,
#' combines them, and lifts over the coordinates from hg19 to hg38 using the UCSC chain file.
#'
#' @param input_file Character vector specifying the path(s) to input file(s).
#'   Each file should be tab-delimited in ChimerDB4 format (or similar).
#' @param output_dir Character string specifying the directory where output files will be written.
#'   Defaults to the current working directory.
#' @param chain_file Character string specifying the path to the chain file. If NULL (default),
#'   the function will download the chain file from UCSC.
#' @param column_names Character vector of column names for the input file(s). If NULL (default),
#'   the function will use the default ChimerDB4 column names.
#'
#' @return A list containing three data frames:
#' \describe{
#'   \item{combined}{Data frame with both hg19 and hg38 coordinates (all files combined)}
#'   \item{hg38_only}{Data frame with only hg38 coordinates for successful liftovers}
#'   \item{failed}{Data frame with records that failed to lift over}
#' }
#'
#' @details The function creates three output files in the specified directory (for *all* files combined):
#' \describe{
#'   \item{chimer_data_hg19_and_hg38.txt}{Contains both hg19 and hg38 coordinates}
#'   \item{chimer_data_hg38.txt}{Contains only hg38 coordinates for successful liftovers}
#'   \item{failed_liftover.txt}{Contains records that failed to lift over}
#' }
#'
#' @examples
#' \dontrun{
#' # Multiple files
#' liftover_fusion_coordinates(
#'   input_file = c("chimerdb4_file1.txt", "chimerdb4_file2.txt"),
#'   output_dir = "path/to/output"
#' )
#'
#' # Single file still works (same interface)
#' liftover_fusion_coordinates("chimerdb4_file.txt")
#' }
#'
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr mutate left_join filter select n
#' @importFrom utils download.file
#' @export
liftover_fusion_coordinates <- function(input_file,
                                        output_dir = ".",
                                        chain_file = NULL,
                                        column_names = NULL) {
  # 1) Possibly process multiple input files
  if (!is.character(input_file)) {
    stop("`input_file` must be a character vector of file paths.")
  }
  if (length(input_file) < 1) {
    stop("`input_file` is empty. Please provide at least one file path.")
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

  # We will read each file, fix columns, and combine into one data frame:
  all_data_list <- list()

  for (f in input_file) {
    if (!file.exists(f)) {
      stop("Input file does not exist: ", f)
    }

    this_data <- read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE)

    # If the file has more columns than our column_names, append column names
    if (ncol(this_data) > length(column_names)) {
      extra_cols <- ncol(this_data) - length(column_names)
      column_names_expanded <- c(
        column_names,
        paste0("V", (length(column_names) + 1):(length(column_names) + extra_cols))
      )
      colnames(this_data) <- column_names_expanded
    } else {
      # If fewer, truncate
      cnames_to_use <- column_names[1:ncol(this_data)]
      colnames(this_data) <- cnames_to_use
    }

    # Store file name (optional) to track where each row came from
    this_data$file_source <- basename(f)

    all_data_list[[f]] <- this_data
  }

  # Combine them all
  chimer_data <- do.call(rbind, all_data_list)

  # Ensure the essential columns exist
  required_cols <- c("Chr1", "Pos1", "Strand1", "Chr2", "Pos2", "Strand2")
  missing_cols <- required_cols[!required_cols %in% colnames(chimer_data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the combined data: ", paste(missing_cols, collapse = ", "))
  }

  # Function to normalize strand
  normalize_strand <- function(strand_values) {
    strand_values <- as.character(strand_values)
    valid_values <- c("+", "-", "*")
    strand_values[!strand_values %in% valid_values] <- "*"
    return(strand_values)
  }

  chimer_data$Strand1 <- normalize_strand(chimer_data$Strand1)
  chimer_data$Strand2 <- normalize_strand(chimer_data$Strand2)

  # 2) Create GRanges objects for each partner
  library(GenomicRanges)
  library(IRanges)
  gr1 <- GRanges(
    seqnames = chimer_data$Chr1,
    ranges   = IRanges(start = as.numeric(chimer_data$Pos1), width = 1),
    strand   = chimer_data$Strand1,
    fusion_id = seq_len(nrow(chimer_data)),
    gene = if ("Gene1" %in% colnames(chimer_data)) chimer_data$Gene1 else NA
  )

  gr2 <- GRanges(
    seqnames = chimer_data$Chr2,
    ranges   = IRanges(start = as.numeric(chimer_data$Pos2), width = 1),
    strand   = chimer_data$Strand2,
    fusion_id = seq_len(nrow(chimer_data)),
    gene = if ("Gene2" %in% colnames(chimer_data)) chimer_data$Gene2 else NA
  )

  # 3) Handle chain file
  if (is.null(chain_file)) {
    chain_file <- file.path(output_dir, "hg19ToHg38.over.chain")
    if (!file.exists(chain_file)) {
      temp_gz <- file.path(output_dir, "hg19ToHg38.over.chain.gz")
      download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                    temp_gz, mode = "wb")
      system(paste("gunzip", shQuote(temp_gz)))
    }
  } else if (!file.exists(chain_file)) {
    stop("Chain file does not exist: ", chain_file)
  }

  library(rtracklayer)
  chain <- import.chain(chain_file)

  # 4) LiftOver
  gr1_hg38 <- liftOver(gr1, chain)
  gr2_hg38 <- liftOver(gr2, chain)

  gr1_hg38_unlisted <- unlist(gr1_hg38)
  gr2_hg38_unlisted <- unlist(gr2_hg38)

  df1_hg38 <- if (length(gr1_hg38_unlisted) > 0) {
    data.frame(
      fusion_id = gr1_hg38_unlisted$fusion_id,
      chr1_hg38 = as.character(seqnames(gr1_hg38_unlisted)),
      pos1_hg38 = start(gr1_hg38_unlisted),
      strand1_hg38 = as.character(strand(gr1_hg38_unlisted)),
      gene1 = gr1_hg38_unlisted$gene
    )
  } else {
    data.frame(
      fusion_id = integer(),
      chr1_hg38 = character(),
      pos1_hg38 = integer(),
      strand1_hg38 = character(),
      gene1 = character()
    )
  }

  df2_hg38 <- if (length(gr2_hg38_unlisted) > 0) {
    data.frame(
      fusion_id = gr2_hg38_unlisted$fusion_id,
      chr2_hg38 = as.character(seqnames(gr2_hg38_unlisted)),
      pos2_hg38 = start(gr2_hg38_unlisted),
      strand2_hg38 = as.character(strand(gr2_hg38_unlisted)),
      gene2 = gr2_hg38_unlisted$gene
    )
  } else {
    data.frame(
      fusion_id = integer(),
      chr2_hg38 = character(),
      pos2_hg38 = integer(),
      strand2_hg38 = character(),
      gene2 = character()
    )
  }

  # 5) Combine lifted results with original table
  library(dplyr)

  # attach df1
  result1 <- chimer_data %>%
    mutate(fusion_id = seq_len(n())) %>%
    left_join(df1_hg38, by = "fusion_id")

  # attach df2
  result2 <- result1 %>%
    left_join(df2_hg38, by = "fusion_id")

  # 6) Fail / success detection
  failed_liftover <- result2 %>%
    filter(is.na(chr1_hg38) | is.na(chr2_hg38))

  if (nrow(failed_liftover) > 0) {
    warning(paste0(nrow(failed_liftover), " record(s) failed to lift over."))
    write.table(
      failed_liftover,
      file.path(output_dir, "failed_liftover.txt"),
      sep="\t", quote=FALSE, row.names=FALSE
    )
  }

  # 7) Final table with both hg19 & hg38
  gene_cols <- intersect(c("Gene1", "Gene2"), colnames(chimer_data))
  base_cols <- setdiff(
    colnames(chimer_data),
    c("Chr1", "Pos1", "Strand1", "Chr2", "Pos2", "Strand2", gene_cols)
  )

  final_result <- result2 %>%
    select(all_of(c(
      base_cols,
      "Gene1", "Chr1", "Pos1", "Strand1", "chr1_hg38", "pos1_hg38", "strand1_hg38",
      "Gene2", "Chr2", "Pos2", "Strand2", "chr2_hg38", "pos2_hg38", "strand2_hg38",
      "file_source"  # so we know which file each row came from
    )))

  write.table(
    final_result,
    file.path(output_dir, "chimer_data_hg19_and_hg38.txt"),
    sep="\t", quote=FALSE, row.names=FALSE
  )

  # 8) hg38-only
  # Construct columns for gene/chr/pos if present
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
  # Add file_source if you like
  select_cols <- c(select_cols, "file_source")

  hg38_only <- result2 %>%
    filter(!is.na(chr1_hg38) & !is.na(chr2_hg38)) %>%
    select(all_of(select_cols))

  # rename columns
  names(hg38_only) <- gsub("_hg38", "", names(hg38_only))
  if ("gene1" %in% names(hg38_only)) names(hg38_only)[names(hg38_only) == "gene1"] <- "Gene1"
  if ("gene2" %in% names(hg38_only)) names(hg38_only)[names(hg38_only) == "gene2"] <- "Gene2"

  write.table(
    hg38_only,
    file.path(output_dir, "chimer_data_hg38.txt"),
    sep="\t", quote=FALSE, row.names=FALSE
  )

  # Summaries
  message("Total records (all files combined): ", nrow(chimer_data))
  message("Successfully lifted over: ", nrow(hg38_only))
  message("Failed to lift over: ", nrow(chimer_data) - nrow(hg38_only))

  return(list(
    combined   = final_result,
    hg38_only  = hg38_only,
    failed     = failed_liftover
  ))
}

#' Lift over fusion gene coordinates from hg19 to hg38 (multi-file version)
#'
#' This function accepts *multiple* input files (a character vector) in ChimerDB4-like format,
#' combines them, and lifts over the coordinates from hg19 to hg38 using the UCSC chain file.
#'
#' @param input_file Character vector specifying the path(s) to input file(s).
#'   Each file should be tab-delimited in ChimerDB4 format (or similar).
#' @param output_dir Character string specifying the directory where output files will be written.
#'   Defaults to the current working directory.
#' @param chain_file Character string specifying the path to the chain file. If NULL (default),
#'   the function will download the chain file from UCSC.
#' @param column_names Character vector of column names for the input file(s). If NULL (default),
#'   the function will use the default ChimerDB4 column names.
#'
#' @return A list containing three data frames:
#' \describe{
#'   \item{combined}{Data frame with both hg19 and hg38 coordinates (all files combined)}
#'   \item{hg38_only}{Data frame with only hg38 coordinates for successful liftovers}
#'   \item{failed}{Data frame with records that failed to lift over}
#' }
#'
#' @details The function creates three output files in the specified directory (for *all* files combined):
#' \describe{
#'   \item{chimer_data_hg19_and_hg38.txt}{Contains both hg19 and hg38 coordinates}
#'   \item{chimer_data_hg38.txt}{Contains only hg38 coordinates for successful liftovers}
#'   \item{failed_liftover.txt}{Contains records that failed to lift over}
#' }
#'
#' @examples
#' \dontrun{
#' # Multiple files
#' liftover_fusion_coordinates(
#'   input_file = c("chimerdb4_file1.txt", "chimerdb4_file2.txt"),
#'   output_dir = "path/to/output"
#' )
#'
#' # Single file still works (same interface)
#' liftover_fusion_coordinates("chimerdb4_file.txt")
#' }
#'
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom dplyr mutate left_join filter select n
#' @importFrom utils download.file
#' @export
liftover_fusion_coordinates <- function(input_file,
                                        output_dir = ".",
                                        chain_file = NULL,
                                        column_names = NULL) {
  # 1) Possibly process multiple input files
  if (!is.character(input_file)) {
    stop("`input_file` must be a character vector of file paths.")
  }
  if (length(input_file) < 1) {
    stop("`input_file` is empty. Please provide at least one file path.")
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

  # We will read each file, fix columns, and combine into one data frame:
  all_data_list <- list()

  for (f in input_file) {
    if (!file.exists(f)) {
      stop("Input file does not exist: ", f)
    }

    this_data <- read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE)

    # If the file has more columns than our column_names, append column names
    if (ncol(this_data) > length(column_names)) {
      extra_cols <- ncol(this_data) - length(column_names)
      column_names_expanded <- c(
        column_names,
        paste0("V", (length(column_names) + 1):(length(column_names) + extra_cols))
      )
      colnames(this_data) <- column_names_expanded
    } else {
      # If fewer, truncate
      cnames_to_use <- column_names[1:ncol(this_data)]
      colnames(this_data) <- cnames_to_use
    }

    # Store file name (optional) to track where each row came from
    this_data$file_source <- basename(f)

    all_data_list[[f]] <- this_data
  }

  # Combine them all
  chimer_data <- do.call(rbind, all_data_list)

  # Ensure the essential columns exist
  required_cols <- c("Chr1", "Pos1", "Strand1", "Chr2", "Pos2", "Strand2")
  missing_cols <- required_cols[!required_cols %in% colnames(chimer_data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns in the combined data: ", paste(missing_cols, collapse = ", "))
  }

  # Function to normalize strand
  normalize_strand <- function(strand_values) {
    strand_values <- as.character(strand_values)
    valid_values <- c("+", "-", "*")
    strand_values[!strand_values %in% valid_values] <- "*"
    return(strand_values)
  }

  chimer_data$Strand1 <- normalize_strand(chimer_data$Strand1)
  chimer_data$Strand2 <- normalize_strand(chimer_data$Strand2)

  # 2) Create GRanges objects for each partner
  library(GenomicRanges)
  library(IRanges)
  gr1 <- GRanges(
    seqnames = chimer_data$Chr1,
    ranges   = IRanges(start = as.numeric(chimer_data$Pos1), width = 1),
    strand   = chimer_data$Strand1,
    fusion_id = seq_len(nrow(chimer_data)),
    gene = if ("Gene1" %in% colnames(chimer_data)) chimer_data$Gene1 else NA
  )

  gr2 <- GRanges(
    seqnames = chimer_data$Chr2,
    ranges   = IRanges(start = as.numeric(chimer_data$Pos2), width = 1),
    strand   = chimer_data$Strand2,
    fusion_id = seq_len(nrow(chimer_data)),
    gene = if ("Gene2" %in% colnames(chimer_data)) chimer_data$Gene2 else NA
  )

  # 3) Handle chain file
  if (is.null(chain_file)) {
    chain_file <- file.path(output_dir, "hg19ToHg38.over.chain")
    if (!file.exists(chain_file)) {
      temp_gz <- file.path(output_dir, "hg19ToHg38.over.chain.gz")
      download.file("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
                    temp_gz, mode = "wb")
      system(paste("gunzip", shQuote(temp_gz)))
    }
  } else if (!file.exists(chain_file)) {
    stop("Chain file does not exist: ", chain_file)
  }

  library(rtracklayer)
  chain <- import.chain(chain_file)

  # 4) LiftOver
  gr1_hg38 <- liftOver(gr1, chain)
  gr2_hg38 <- liftOver(gr2, chain)

  gr1_hg38_unlisted <- unlist(gr1_hg38)
  gr2_hg38_unlisted <- unlist(gr2_hg38)

  df1_hg38 <- if (length(gr1_hg38_unlisted) > 0) {
    data.frame(
      fusion_id = gr1_hg38_unlisted$fusion_id,
      chr1_hg38 = as.character(seqnames(gr1_hg38_unlisted)),
      pos1_hg38 = start(gr1_hg38_unlisted),
      strand1_hg38 = as.character(strand(gr1_hg38_unlisted)),
      gene1 = gr1_hg38_unlisted$gene
    )
  } else {
    data.frame(
      fusion_id = integer(),
      chr1_hg38 = character(),
      pos1_hg38 = integer(),
      strand1_hg38 = character(),
      gene1 = character()
    )
  }

  df2_hg38 <- if (length(gr2_hg38_unlisted) > 0) {
    data.frame(
      fusion_id = gr2_hg38_unlisted$fusion_id,
      chr2_hg38 = as.character(seqnames(gr2_hg38_unlisted)),
      pos2_hg38 = start(gr2_hg38_unlisted),
      strand2_hg38 = as.character(strand(gr2_hg38_unlisted)),
      gene2 = gr2_hg38_unlisted$gene
    )
  } else {
    data.frame(
      fusion_id = integer(),
      chr2_hg38 = character(),
      pos2_hg38 = integer(),
      strand2_hg38 = character(),
      gene2 = character()
    )
  }

  # 5) Combine lifted results with original table
  library(dplyr)

  # attach df1
  result1 <- chimer_data %>%
    mutate(fusion_id = seq_len(n())) %>%
    left_join(df1_hg38, by = "fusion_id")

  # attach df2
  result2 <- result1 %>%
    left_join(df2_hg38, by = "fusion_id")

  # 6) Fail / success detection
  failed_liftover <- result2 %>%
    filter(is.na(chr1_hg38) | is.na(chr2_hg38))

  if (nrow(failed_liftover) > 0) {
    warning(paste0(nrow(failed_liftover), " record(s) failed to lift over."))
    write.table(
      failed_liftover,
      file.path(output_dir, "failed_liftover.txt"),
      sep="\t", quote=FALSE, row.names=FALSE
    )
  }

  # 7) Final table with both hg19 & hg38
  gene_cols <- intersect(c("Gene1", "Gene2"), colnames(chimer_data))
  base_cols <- setdiff(
    colnames(chimer_data),
    c("Chr1", "Pos1", "Strand1", "Chr2", "Pos2", "Strand2", gene_cols)
  )

  final_result <- result2 %>%
    select(all_of(c(
      base_cols,
      "Gene1", "Chr1", "Pos1", "Strand1", "chr1_hg38", "pos1_hg38", "strand1_hg38",
      "Gene2", "Chr2", "Pos2", "Strand2", "chr2_hg38", "pos2_hg38", "strand2_hg38",
      "file_source"  # so we know which file each row came from
    )))

  write.table(
    final_result,
    file.path(output_dir, "chimer_data_hg19_and_hg38.txt"),
    sep="\t", quote=FALSE, row.names=FALSE
  )

  # 8) hg38-only
  # Construct columns for gene/chr/pos if present
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
  # Add file_source if you like
  select_cols <- c(select_cols, "file_source")

  hg38_only <- result2 %>%
    filter(!is.na(chr1_hg38) & !is.na(chr2_hg38)) %>%
    select(all_of(select_cols))

  # rename columns
  names(hg38_only) <- gsub("_hg38", "", names(hg38_only))
  if ("gene1" %in% names(hg38_only)) names(hg38_only)[names(hg38_only) == "gene1"] <- "Gene1"
  if ("gene2" %in% names(hg38_only)) names(hg38_only)[names(hg38_only) == "gene2"] <- "Gene2"

  write.table(
    hg38_only,
    file.path(output_dir, "chimer_data_hg38.txt"),
    sep="\t", quote=FALSE, row.names=FALSE
  )

  # Summaries
  message("Total records (all files combined): ", nrow(chimer_data))
  message("Successfully lifted over: ", nrow(hg38_only))
  message("Failed to lift over: ", nrow(chimer_data) - nrow(hg38_only))

  return(list(
    combined   = final_result,
    hg38_only  = hg38_only,
    failed     = failed_liftover
  ))
}

#' Lift over a vector of genomic coordinates
#'
#' Given a character vector of coordinates in the format "chrN:start-end", this function
#' uses the UCSC \code{liftOver} utility (via the \code{rtracklayer} package) to map them
#' to another assembly, for example from hg19 to hg38.
#'
#' @param coords A character vector of coordinate strings like \code{c("chr1:100-400", "chr3:450000-9890098")}.
#' @param chain_file Path to the chain file for liftOver (e.g. "hg19ToHg38.over.chain"). If this
#'   file does not exist and \code{download_chain=TRUE}, it will be downloaded.
#' @param download_chain Logical, if \code{TRUE} and \code{chain_file} does not exist,
#'   attempts to download the hg19->hg38 chain file from UCSC to the specified path.
#' @param genome_style Logical, if \code{TRUE}, tries to maintain "chr" prefix. If your chain
#'   file expects a certain style, ensure your input \code{coords} have matching "chr" or no-"chr"
#'   style. Alternatively, you can remove or add the prefix before calling this function.
#' @return A data frame with the following columns:
#'   \describe{
#'     \item{original}{The original coordinate string.}
#'     \item{mapped_chr}{The chromosome (in new assembly).}
#'     \item{mapped_start}{The start position in new assembly.}
#'     \item{mapped_end}{The end position in new assembly.}
#'     \item{status}{Either "mapped" if a single unique mapping was found, "failed" if no mapping
#'       was found, or "split" if the range mapped to multiple segments.}
#'   }
#'
#' @examples
#' \dontrun{
#' liftover_coordinates(c("chr22:23632600-23632600", "chr9:133729451-133729600"),
#'                      chain_file="hg19ToHg38.over.chain")
#' }
#'
#' @importFrom utils download.file
#' @importFrom rtracklayer import.chain liftOver
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @export
liftover_coordinates <- function(coords,
                                 chain_file = "hg19ToHg38.over.chain",
                                 download_chain = FALSE,
                                 genome_style = TRUE) {
  if (!requireNamespace("rtracklayer", quietly=TRUE))
    stop("Package 'rtracklayer' is required but not installed.")

  if (!requireNamespace("GenomicRanges", quietly=TRUE))
    stop("Package 'GenomicRanges' is required but not installed.")

  if (!requireNamespace("IRanges", quietly=TRUE))
    stop("Package 'IRanges' is required but not installed.")

  # If chain file doesn't exist, optionally download from UCSC
  if (!file.exists(chain_file)) {
    if (download_chain) {
      message("Downloading chain file from UCSC...")
      tmp_gz <- paste0(chain_file, ".gz")
      utils::download.file(
        url = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz",
        destfile = tmp_gz,
        mode = "wb"
      )
      system(paste("gunzip -c", shQuote(tmp_gz), ">", shQuote(chain_file)))
      unlink(tmp_gz)
      if (!file.exists(chain_file)) {
        stop("Failed to download or decompress chain file.")
      }
    } else {
      stop("Chain file does not exist: ", chain_file)
    }
  }

  # Parse the input coords into a data frame
  parse_coord <- function(x) {
    # "chr3:123-456"
    # contig = chr3
    # range_part = 123-456
    contig <- sub(":.*", "", x)
    range_part <- sub(".*:", "", x)
    start_str <- sub("-.*", "", range_part)
    end_str   <- sub(".*-", "", range_part)
    start_num <- as.numeric(start_str)
    end_num   <- as.numeric(end_str)

    # If user wrote e.g. "chr1:100" with no dash, that will produce NA for end.
    # We can handle that by setting end = start if so.
    if (is.na(end_num)) {
      end_num <- start_num
    }

    list(
      original = x,
      seqnames = contig,
      start    = min(start_num, end_num),
      end      = max(start_num, end_num)
    )
  }

  parsed <- lapply(coords, parse_coord)
  df_parsed <- do.call(rbind, lapply(parsed, as.data.frame))
  rownames(df_parsed) <- NULL

  # If user wants to ensure "chr" prefix, do so here:
  if (genome_style) {
    # We'll assume that if it doesn't start with "chr", add it:
    df_parsed$seqnames <- ifelse(grepl("^chr", df_parsed$seqnames, ignore.case=TRUE),
                                 df_parsed$seqnames,
                                 paste0("chr", df_parsed$seqnames))
  }

  # Create GRanges
  library(GenomicRanges)
  library(IRanges)
  gr_input <- GRanges(
    seqnames = df_parsed$seqnames,
    ranges   = IRanges(start=df_parsed$start, end=df_parsed$end),
    original = df_parsed$original
  )

  # Import chain, run liftOver
  library(rtracklayer)
  chain <- import.chain(chain_file)

  lifted <- liftOver(gr_input, chain)

  # Build a results data frame
  # Each element in `lifted` is a GRanges of length 0 or >1 if multi-mapping
  output <- data.frame(
    original     = df_parsed$original,
    mapped_chr   = NA_character_,
    mapped_start = NA_integer_,
    mapped_end   = NA_integer_,
    status       = "failed",
    stringsAsFactors = FALSE
  )

  for (i in seq_along(lifted)) {
    if (length(lifted[[i]]) == 0) {
      # no mapping
      output$status[i] <- "failed"
    } else if (length(lifted[[i]]) == 1) {
      # single mapping
      output$mapped_chr[i]   <- as.character(seqnames(lifted[[i]]))[1]
      output$mapped_start[i] <- start(lifted[[i]])[1]
      output$mapped_end[i]   <- end(lifted[[i]])[1]
      output$status[i]       <- "mapped"
    } else {
      # multiple intervals
      # (for example, if a large region spanned a gap)
      # We'll record the first one, or you can do something else
      output$mapped_chr[i]   <- as.character(seqnames(lifted[[i]]))[1]
      output$mapped_start[i] <- start(lifted[[i]])[1]
      output$mapped_end[i]   <- end(lifted[[i]])[1]
      output$status[i]       <- "split"
      # Or you could keep them all if you prefer a more advanced approach
    }
  }

  return(output)
}


