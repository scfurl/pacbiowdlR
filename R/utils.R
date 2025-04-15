#' Annotate Coverage Data with Gene Information
#'
#' Converts a coverage data frame into a GRanges object, extracts gene ranges from a transcript
#' database, finds overlaps between the coverage intervals and genes, and maps gene IDs to gene
#' symbols (or returns Entrez IDs).
#'
#' @param coverage_df A data frame with at least the following columns: \code{seqnames},
#'   \code{start}, and \code{end} representing coverage intervals.
#' @param txdb A transcript database object from which gene information is extracted.
#'   Defaults to \code{TxDb.Hsapiens.UCSC.hg19.knownGene}.
#' @param return_entrez Logical. If \code{FALSE} (default), gene symbols are returned; otherwise,
#'   the Entrez gene IDs are returned.
#'
#' @return The original \code{coverage_df} with two additional columns: \code{gene_id} (gene IDs)
#'   and \code{gene_symbol} (gene symbols or Entrez IDs, depending on \code{return_entrez}).
#'
#' @details The function converts the input data frame to a \code{GRanges} object using
#'   \code{GRanges()} and \code{IRanges()}. It then obtains gene ranges via \code{genes()} from the
#'   transcript database and finds overlaps with \code{findOverlaps()}. Gene IDs are extracted from the
#'   metadata columns of the gene ranges, and, if requested, mapped to gene symbols using
#'   \code{mapIds()}.
#'
#' @importFrom GenomicRanges GRanges findOverlaps
#' @importFrom GenomicFeatures genes
#' @importFrom IRanges IRanges
#' @importFrom dplyr mutate lag
#' @importFrom AnnotationDbi mapIds
#' @keywords internal
#' @export
annotateCoverageWithGenes <- function(coverage_df,
                                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      return_entrez = FALSE) {
  # Convert coverage data frame into GRanges.
  coverage_gr <- GRanges(
    seqnames = coverage_df$seqnames,
    ranges   = IRanges(start = coverage_df$start, end = coverage_df$end)
  )

  # Extract gene ranges from the transcript database.
  gene_ranges <- genes(txdb, single.strand.genes.only = FALSE)

  # Find overlaps between coverage intervals and gene ranges.
  ov <- findOverlaps(coverage_gr, gene_ranges)

  # Map to gene IDs.
  coverage_df$gene_id <- NA_character_
  coverage_df$gene_id[queryHits(ov)] <- mcols(gene_ranges)$gene_id[subjectHits(ov)]

  if (!return_entrez) {
    # Convert Entrez IDs to gene symbols.
    suppressMessages({
      coverage_df$gene_symbol <- mapIds(
        org.Hs.eg.db,
        keys      = gsub("^GeneID:", "", coverage_df$gene_id),  # Remove "GeneID:" prefix, if present.
        column    = "SYMBOL",
        keytype   = "ENTREZID",
        multiVals = "first"
      )
    })
  } else {
    coverage_df$gene_symbol <- coverage_df$gene_id
  }

  coverage_df
}

#' Generate Cumulative Genomic Offsets
#'
#' Computes cumulative genomic positions for coverage data and variant calls by calculating
#' chromosome offsets.
#'
#' @param coverage_data A data frame with columns \code{seqnames}, \code{start}, and \code{end}
#'   representing coverage intervals.
#' @param variant_calls A data frame with columns \code{seqnames}, \code{start}, and \code{end}
#'   representing variant calls.
#'
#' @return A list with two elements:
#'   \describe{
#'     \item{\code{coverage_data}}{The input coverage data augmented with a column
#'       \code{cumulative_genomic_position}.}
#'     \item{\code{variant_calls}}{The input variant calls augmented with columns
#'       \code{cumulative_start} and \code{cumulative_end}.}
#'   }
#'
#' @details The function filters for valid chromosomes (chr1 - chr22, X, Y), calculates the maximum
#'   \code{end} of each chromosome as a proxy for its length, and uses the cumulative sum of chromosome
#'   lengths to determine an offset for each chromosome. It then adds this offset to the \code{start} of
#'   each interval to obtain a cumulative genomic position.
#'
#' @importFrom dplyr filter mutate left_join group_by summarize arrange
#' @keywords internal
#' @export
generate_offsets <- function(coverage_data, variant_calls) {
  # Define the chromosome order and valid chromosomes.
  chromosome_order <- c(as.character(1:22), "X", "Y")
  valid_chroms <- paste0("chr", chromosome_order)

  # Filter and factorize seqnames.
  coverage_data <- coverage_data %>%
    filter(seqnames %in% valid_chroms) %>%
    mutate(seqnames = factor(seqnames, levels = valid_chroms))

  # Compute chromosome lengths (using max(end) as a proxy) and offsets.
  chr_offsets <- coverage_data %>%
    group_by(seqnames) %>%
    summarize(chr_length = max(end), .groups = "drop") %>%  # using max(end) as a proxy for length
    arrange(factor(seqnames, levels = valid_chroms)) %>%
    mutate(offset = lag(cumsum(as.numeric(chr_length)), default = 0))

  # Add cumulative positions for coverage data.
  coverage_data <- coverage_data %>%
    left_join(chr_offsets, by = "seqnames") %>%
    mutate(cumulative_genomic_position = as.numeric(start) + offset)

  # Add cumulative positions for variant calls.
  variant_calls <- variant_calls %>%
    filter(seqnames %in% valid_chroms) %>%
    left_join(chr_offsets, by = "seqnames") %>%
    mutate(cumulative_start = as.numeric(start) + offset,
           cumulative_end   = as.numeric(end) + offset)

  return(list(coverage_data = coverage_data, variant_calls = variant_calls))
}

#' Generate CNA Plot from Regions
#'
#' Generates a genome-wide copy number alteration (CNA) plot by combining coverage data from a BigWig
#' file and variant calls from a VCF file. Delta values are computed from the coverage data and overlaid
#' with variant regions. Optionally, cytoband annotations can be added.
#'
#' @param depth_bigwig_file Character. File path to a BigWig file containing coverage (depth) information.
#' @param variant_file Character. File path to a VCF file containing variant call information.
#' @param threshold Numeric. Delta threshold used for display. Default is 1.5.
#' @param downsample Numeric. Proportion of the coverage data to retain after downsampling. Default is 0.01.
#' @param point_size Numeric. Size of the points in the plot. Default is 0.5.
#' @param colors Named vector. Colors to be used for chromosomes; if \code{NULL}, default colors are assigned.
#' @param log_scale Logical. If \code{TRUE}, applies a log scale (unused in code, but available for future extension). Default is \code{TRUE}.
#' @param remove_centromeres Logical. If \code{TRUE}, centromere regions are removed (unused in code, but available for extension). Default is \code{TRUE}.
#' @param max_value Numeric. Maximum delta value for setting plot y-axis limits. Default is \code{NULL}.
#' @param min_value Numeric. Minimum delta value for setting plot y-axis limits. Default is \code{NULL}.
#' @param min_variant_distance Numeric. Minimum distance (in bp) for a variant call to be retained. Default is 10000.
#' @param samplename Character. A label for the sample which is displayed in the plot title.
#' @param show_cytobands Logical. If \code{TRUE}, cytoband labels are added to the variant calls. Default is \code{FALSE}.
#' @param cytoband_angle Numeric. Angle (in degrees) to rotate the cytoband text labels. Default is 90.
#' @param chr_filter Character. If provided, limits the analysis to the specified chromosome.
#' @param return_data Logical. If \code{TRUE}, returns a list containing the plot along with the underlying
#'   coverage and variant call data. Default is \code{FALSE}.
#'
#' @return Either a \code{ggplot2} object representing the CNA plot, or, if \code{return_data} is \code{TRUE},
#'   a list with elements \code{plot}, \code{coverage_data}, and \code{variant_calls}.
#'
#' @details This function performs the following steps:
#'   \enumerate{
#'     \item Imports coverage data from a BigWig file using \code{import.bw()} from \code{rtracklayer}
#'           and converts it into a data frame.
#'     \item Filters for valid chromosomes (chr1-22, X, Y), optionally restricting to a single chromosome.
#'     \item Downsamples the coverage data and computes a delta value for each data point.
#'     \item Reads variant calls from a VCF file using \code{fread()} from \code{data.table}, parsing out
#'           the end positions and variant types.
#'     \item Optionally annotates cytoband information for each variant call using overlaps with
#'           \code{hg19IdeogramCyto} (assumed to be available in the global environment).
#'     \item Filters variant calls by a minimum distance if specified.
#'     \item Computes cumulative genomic offsets via the \code{generate_offsets()} function.
#'     \item Constructs and returns a plot built with \code{ggplot2} that shows the coverage data points,
#'           highlighted variant regions, and optionally cytoband labels.
#'   }
#'
#' @importFrom rtracklayer import.bw
#' @importFrom dplyr filter mutate left_join group_by summarize arrange
#' @importFrom data.table fread
#' @importFrom ggplot2 ggplot geom_point geom_rect labs theme_minimal theme guides scale_color_manual scale_fill_manual scale_y_continuous geom_text
#' @keywords internal
#' @export
generateCNAPlotFromRegions2 <- function(depth_bigwig_file, variant_file, threshold = 1.5,
                                        downsample = 0.01, point_size = 0.5, colors = NULL,
                                        log_scale = TRUE, remove_centromeres = TRUE,
                                        max_value = NULL, min_value = NULL,
                                        min_variant_distance = 10000, samplename = "",
                                        show_cytobands = FALSE, cytoband_angle = 90,
                                        chr_filter = NULL, return_data = FALSE) {
  # Import coverage data from BigWig file.
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- as.data.frame(coverage_data)

  # Filter for valid chromosomes.
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))
  coverage_data <- coverage_data %>% filter(seqnames %in% valid_chromosomes)

  # Optionally, filter for a specific chromosome.
  if (!is.null(chr_filter)) {
    coverage_data <- coverage_data %>% filter(seqnames == chr_filter)
  }

  # Downsample the coverage data.
  set.seed(42)
  sampled_indices <- sample(seq_len(nrow(coverage_data)),
                            size = max(1, floor(nrow(coverage_data) * downsample)))
  coverage_data <- coverage_data[sampled_indices, ]

  # Compute the delta (difference from the mean coverage).
  mean_coverage <- mean(coverage_data$score)
  coverage_data$delta <- coverage_data$score - mean_coverage

  # Read variant calls from the VCF file.
  vcf <- fread(variant_file, skip = "#CHROM")
  colnames(vcf)[1] <- "CHROM"
  vcf <- vcf %>% filter(CHROM %in% valid_chromosomes)

  # Optionally, filter variant calls by a specific chromosome.
  if (!is.null(chr_filter)) {
    vcf <- vcf %>% filter(CHROM == chr_filter)
  }

  # Extract variant information.
  end <- strsplit(vcf$INFO, ";") %>% sapply("[[", 3) %>% gsub("END=", "", .) %>% as.numeric()
  type <- strsplit(vcf$INFO, ";") %>% sapply("[[", 2) %>% gsub("SVTYPE=", "", .)
  variant_calls <- data.frame(seqnames = vcf$CHROM, start = vcf$POS, end = end, type = type)

  # Convert variant_calls to a GRanges object.
  variant_gr <- GRanges(
    seqnames = variant_calls$seqnames,
    ranges = IRanges(start = variant_calls$start, end = variant_calls$end)
  )

  # Optionally, annotate cytoband information.
  if (show_cytobands) {
    overlaps <- findOverlaps(variant_gr, hg19IdeogramCyto, type = "any")
    query_indices <- queryHits(overlaps)
    subject_indices <- subjectHits(overlaps)
    overlap_df <- data.frame(
      query = query_indices,
      cytoband = hg19IdeogramCyto$name[subject_indices]
    )

    cytoband_annotation <- overlap_df %>%
      mutate(cytoband = as.character(cytoband)) %>%
      group_by(query) %>%
      summarize(
        cytoband_name = if (n() == 1) dplyr::first(cytoband) else paste(dplyr::first(cytoband), dplyr::last(cytoband), sep = "-"),
        .groups = "drop"
      )

    variant_calls$cytoband <- cytoband_annotation$cytoband_name[match(1:nrow(variant_calls), cytoband_annotation$query)]
  } else {
    variant_calls$cytoband <- NA
  }

  # Filter variant calls by minimum distance.
  if (!is.null(min_variant_distance) && min_variant_distance > 0) {
    variant_calls <- variant_calls %>%
      mutate(distance = end - start) %>%
      filter(distance >= min_variant_distance)
  }

  # Generate cumulative offsets.
  offset <- generate_offsets(coverage_data = coverage_data, variant_calls = variant_calls)

  # Set default colors if not provided.
  if (is.null(colors)) {
    chromosome_order <- c(as.character(1:22), "X", "Y")
    colors <- rep(c("black", "gray"), length.out = length(chromosome_order))
  }

  # Create data for variant highlight rectangles.
  if (nrow(offset$variant_calls) > 0) {
    data_highlight <- data.frame(
      xmin = offset$variant_calls$cumulative_start,
      xmax = offset$variant_calls$cumulative_end,
      ymin = rep(-Inf, nrow(offset$variant_calls)),
      ymax = rep(Inf, nrow(offset$variant_calls)),
      VariantType = offset$variant_calls$type,
      cytoband = offset$variant_calls$cytoband,
      seqnames = offset$variant_calls$seqnames,
      x_center = (offset$variant_calls$cumulative_start + offset$variant_calls$cumulative_end) / 2
    )
    data_highlight$label <- paste(data_highlight$seqnames, data_highlight$cytoband, sep = ": ")
    data_highlight$text_color <- ifelse(data_highlight$VariantType == "DEL", "red",
                                        ifelse(data_highlight$VariantType == "DUP", "green", "black"))
  } else {
    data_highlight <- data.frame(
      xmin = numeric(0),
      xmax = numeric(0),
      ymin = numeric(0),
      ymax = numeric(0),
      VariantType = character(0),
      cytoband = character(0),
      seqnames = character(0),
      x_center = numeric(0),
      label = character(0),
      text_color = character(0)
    )
  }

  # Build the plot.
  p <- ggplot() +
    geom_point(data = offset$coverage_data,
               aes(x = cumulative_genomic_position, y = delta, color = seqnames),
               size = point_size) +
    scale_color_manual(values = colors) +
    geom_rect(data = data_highlight,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = VariantType),
              color = "transparent",
              alpha = 0.2) +
    scale_fill_manual(values = c("DEL" = "red", "DUP" = "green")) +
    labs(title = paste0("Genome-Wide CNA ", samplename),
         x = NULL,
         y = "Delta") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    guides(color = "none") +
    scale_y_continuous(limits = c(
      ifelse(is.null(min_value), min(offset$coverage_data$delta, na.rm = TRUE), min_value),
      ifelse(is.null(max_value), max(offset$coverage_data$delta, na.rm = TRUE), max_value)
    ))

  # Optionally add cytoband labels.
  if (show_cytobands) {
    label_y <- if (!is.null(max_value)) {
      max_value * 0.9
    } else {
      max(coverage_data$delta, na.rm = TRUE) * 0.9
    }
    data_highlight_labels <- data_highlight %>% filter(!is.na(cytoband))
    p <- p + geom_text(
      data = data_highlight_labels,
      aes(x = x_center, y = label_y, label = label, color = I(text_color)),
      size = 3, angle = cytoband_angle, hjust = 0
    )
  }

  # Return the plot (and optionally, additional data).
  if (return_data) {
    return(list(plot = p,
                coverage_data = offset$coverage_data,
                variant_calls = offset$variant_calls))
  } else {
    return(p)
  }
}
