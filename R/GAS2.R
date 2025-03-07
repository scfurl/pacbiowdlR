#' @export
createGeneRegions <- function(
    genes = NULL,
    useTSS = FALSE,
    extendTSS = TRUE,
    geneUpstream = 5000,
    geneDownstream = 0,
    geneScaleFactor = 5,
    includeChr = NULL  # Changed from excludeChr to includeChr
) {

  # Validate inputs
  if(is.null(genes)) {
    stop("Please provide a GRanges object containing gene annotations!")
  }

  if(!inherits(genes, "GRanges")) {
    stop("genes must be a GRanges object!")
  }

  if(inherits(mcols(genes)$symbol, "list") | inherits(mcols(genes)$symbol, "SimpleList")) {
    stop("Found a list in genes symbol! This is an incorrect format. Please correct your genes!")
  }

  if(!any(colnames(mcols(genes)) == "symbol")) {
    stop("No symbol column in genes! A column named symbol is expected in the GRanges object!")
  }

  # Filter genes by included chromosomes if specified
  if(!is.null(includeChr)) {
    geneRegions <- genes[as.character(seqnames(genes)) %in% includeChr,]
    seqlevels(geneRegions) <- as.character(unique(seqnames(geneRegions)))
  } else {
    # If no includeChr is specified, use all chromosomes
    geneRegions <- genes
  }

  # Remove genes with NA symbols
  geneRegions <- geneRegions[!is.na(mcols(geneRegions)$symbol)]

  # Create Gene Regions based on method
  if(useTSS) {
    message("Creating gene regions using TSS as reference")
    distMethod <- "GenePromoter"

    # Store the original gene start and end positions
    geneRegions$geneStart <- start(GenomicRanges::resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(GenomicRanges::resize(geneRegions, 1, "end"))

    # Resize to just the TSS
    geneRegions <- GenomicRanges::resize(geneRegions, 1, "start")

    # Extend TSS if requested
    if(extendTSS) {
      geneRegions <- extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream)
    }

    # Assign gene weights (all equal for TSS-based method)
    geneRegions$geneWeight <- geneScaleFactor

  } else {
    message("Creating gene regions using gene body as reference")
    distMethod <- "GeneBody"

    # Store the original gene start and end positions
    geneRegions$geneStart <- start(GenomicRanges::resize(geneRegions, 1, "start"))
    geneRegions$geneEnd <- start(GenomicRanges::resize(geneRegions, 1, "end"))

    # Extend gene body
    geneRegions <- extendGR(gr = geneRegions, upstream = geneUpstream, downstream = geneDownstream)

    # Create scaling factor based on gene length
    # Shorter genes will have higher weights
    m <- 1 / width(geneRegions)
    geneRegions$geneWeight <- 1 + m * (geneScaleFactor - 1) / (max(m) - min(m))
  }

  # Sort the gene regions for consistency
  geneRegions <- sort(sortSeqlevels(geneRegions), ignore.strand = TRUE)

  # # Split by chromosome and add index
  # geneRegionsByChr <- split(geneRegions, seqnames(geneRegions))
  # geneRegionsByChr <- lapply(geneRegionsByChr, function(x) {
  #   mcols(x)$idx <- seq_along(x)
  #   return(x)
  # })

  # return(list(
  #   geneRegions = geneRegions,           # Full set of gene regions
  #   geneRegionsByChr = geneRegionsByChr, # Split by chromosome with indices
  #   distMethod = distMethod              # Method used for distance calculation
  # ))
  geneRegions
}

# Helper function to extend GRanges objects
extendGR <- function(gr, upstream = 0, downstream = 0) {
  if(upstream == 0 && downstream == 0) {
    return(gr)
  }

  strand_is_minus <- strand(gr) == "-"

  # For minus strand, upstream and downstream are reversed
  start_shift <- ifelse(strand_is_minus, -downstream, -upstream)
  end_shift <- ifelse(strand_is_minus, upstream, downstream)

  # Adjust start and end positions based on strand
  new_start <- start(gr) + start_shift
  new_end <- end(gr) + end_shift

  # Ensure start positions are not less than 1
  new_start <- pmax(1, new_start)

  # Create new GRanges object
  new_gr <- GRanges(
    seqnames = seqnames(gr),
    ranges = IRanges(start = new_start, end = new_end),
    strand = strand(gr),
    mcols(gr)
  )

  return(new_gr)
}

# Helper function to sort seqlevels in a consistent manner
sortSeqlevels <- function(gr) {
  current_levels <- seqlevels(gr)

  # Common chromosome naming patterns
  std_chromosomes <- paste0("chr", c(1:22, "X", "Y", "M"))

  # Try to sort chromosomes in a standard way
  sorted_levels <- current_levels[order(match(current_levels, std_chromosomes, nomatch = length(std_chromosomes) + 1))]

  # Apply the sorted levels
  seqlevels(gr) <- sorted_levels
  return(gr)
}

#' Process FIRE Peak Scores Over Pre-Mapped Gene Regions
#'
#' This function takes FIRE peak scores already mapped to gene regions
#' and applies a weighting approach similar to ArchR's gene score calculation.
#'
#' @param mappedBedFile Path to the BED file with FIRE peaks mapped to gene regions
#' @param outputFile Path to save the processed scores (default: NULL)
#' @param geneModel A string giving a model function used for weighting peaks. This string
#'   should be a function of `x`, where `x` is the distance from the reference point
#' @param distanceColumn Name of the column containing the distance information (default: "distance")
#' @param scoreColumn Name of the column containing the original FIRE scores (default: "score")
#' @param geneIdColumn Name of the column containing gene identifiers (default: "name")
#' @param geneLengthColumn Name of the column containing gene lengths (default: NULL)
#' @param geneScaleFactor Numeric scaling factor to weight genes based on inverse length
#' @param useTSS Boolean indicating whether distances are relative to TSS (TRUE) or gene body (FALSE)
#'
#' @return A data frame with processed FIRE scores
#' @export
processMappedFIREScores <- function(
    mappedBedFile,
    outputFile = NULL,
    geneModel = "exp(-abs(x)/5000) + exp(-1)",
    distanceColumn = "distance",
    scoreColumn = "score",
    geneIdColumn = "name",
    geneLengthColumn = NULL,
    geneScaleFactor = 5,
    useTSS = FALSE
) {
  # Load required libraries
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
  })

  # Read the mapped BED file
  message("Reading pre-mapped FIRE scores...")
  if(is.character(mappedBedFile)) {
    # Read file if string path provided
    if(grepl("\\.bed$", mappedBedFile)) {
      # Read as BED format
      peaks <- read.table(mappedBedFile, header = FALSE, stringsAsFactors = FALSE)
      # Assign column names based on BED format
      if(ncol(peaks) >= 6) {
        colnames(peaks)[1:6] <- c("chrom", "start", "end", "name", "score", "strand")
      } else if(ncol(peaks) >= 5) {
        colnames(peaks)[1:5] <- c("chrom", "start", "end", "name", "score")
      } else if(ncol(peaks) >= 3) {
        colnames(peaks)[1:3] <- c("chrom", "start", "end")
      }

      # If we have extra columns, try to find the distance column
      if(ncol(peaks) > 6) {
        # Check if distance column is already present
        if(!distanceColumn %in% colnames(peaks)) {
          # Try to find a column that might contain distance information
          # This is just a guess - user might need to specify the correct column
          message("Distance column not explicitly found, attempting to infer...")
          if(all(is.numeric(peaks[[7]]))) {
            colnames(peaks)[7] <- distanceColumn
          } else {
            stop("Could not determine distance column. Please specify the correct column name.")
          }
        }
      } else {
        stop("The BED file does not contain enough columns to process distances.")
      }

    } else if(grepl("\\.csv$", mappedBedFile)) {
      # Read as CSV
      peaks <- read.csv(mappedBedFile, stringsAsFactors = FALSE)
    } else if(grepl("\\.tsv$|\\.txt$", mappedBedFile)) {
      # Read as TSV
      peaks <- read.table(mappedBedFile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    } else {
      stop("Unrecognized file format. Please use .bed, .csv, or .tsv/.txt")
    }
  } else if(is.data.frame(mappedBedFile)) {
    # Use the data frame directly
    peaks <- mappedBedFile
  } else {
    stop("mappedBedFile must be a file path or data frame")
  }

  # Verify required columns exist
  required_cols <- c(distanceColumn, scoreColumn, geneIdColumn)
  missing_cols <- required_cols[!required_cols %in% colnames(peaks)]

  if(length(missing_cols) > 0) {
    stop(paste("Missing required columns:", paste(missing_cols, collapse = ", ")))
  }

  # Convert to data.table for faster processing
  peaks <- as.data.table(peaks)

  # Apply gene weights based on gene length if available
  if(!is.null(geneLengthColumn) && geneLengthColumn %in% colnames(peaks)) {
    message("Calculating gene weights based on gene length...")
    if(useTSS) {
      # For TSS-based approach, all genes get the same weight
      peaks$geneWeight <- geneScaleFactor
    } else {
      # For gene body approach, weight inversely proportional to gene length
      lengths <- peaks[[geneLengthColumn]]

      # Calculate inverse proportional weight
      if(length(unique(lengths)) > 1) {
        m <- 1 / lengths
        peaks$geneWeight <- 1 + (m - min(m)) * (geneScaleFactor - 1) / (max(m) - min(m))
      } else {
        # All genes have same length
        peaks$geneWeight <- geneScaleFactor
      }
    }
  } else {
    # If no gene length column, use a default weight
    message("No gene length column provided, using default gene weight...")
    peaks$geneWeight <- geneScaleFactor
  }

  # Apply gene model to calculate weights based on distance
  message("Applying gene model to calculate weights...")

  # Set x to be the distances for use in the gene model formula
  x <- peaks[[distanceColumn]]

  # Evaluate the gene model expression
  peaks$model_weight <- eval(parse(text = geneModel))

  # Apply gene weights to model weights
  peaks$final_weight <- peaks$model_weight * peaks$geneWeight

  # Calculate weighted scores
  message("Calculating weighted scores...")
  peaks$weighted_score <- peaks[[scoreColumn]] * peaks$final_weight

  # Summarize by gene if needed
  if(length(unique(peaks[[geneIdColumn]])) < nrow(peaks)) {
    message("Multiple peaks per gene detected, summarizing scores by gene...")
    gene_summary <- peaks[, .(
      total_score = sum(get(scoreColumn)),
      weighted_score = sum(weighted_score),
      peak_count = .N,
      max_weight = max(final_weight),
      min_weight = min(final_weight),
      avg_weight = mean(final_weight)
    ), by = geneIdColumn]

    # Merge the summary back with original data if desired
    peaks <- merge(peaks, gene_summary, by = geneIdColumn)
  }

  # Save results if outputFile is provided
  if(!is.null(outputFile)) {
    message(paste0("Saving results to: ", outputFile))

    if(grepl("\\.csv$", outputFile)) {
      write.csv(peaks, outputFile, row.names = FALSE)
    } else if(grepl("\\.tsv$|\\.txt$", outputFile)) {
      write.table(peaks, outputFile, row.names = FALSE, sep = "\t", quote = FALSE)
    } else {
      write.table(peaks, outputFile, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }

  return(as.data.frame(peaks))
}
