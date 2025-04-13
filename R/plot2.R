
library(ggplot2)
library(rtracklayer)
library(dplyr)
library(GenomicRanges)
library(data.table)

#' @export
annotateCoverageWithGenes <- function(coverage_df,
                                      txdb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                                      return_entrez = FALSE) {
  # Convert coverage data frame into GRanges
  coverage_gr <- GRanges(
    seqnames = coverage_df$seqnames,
    ranges   = IRanges(start = coverage_df$start, end = coverage_df$end)
  )

  # Extract gene ranges from TxDb
  gene_ranges <- genes(txdb, single.strand.genes.only=FALSE)

  # Find overlaps
  ov <- findOverlaps(coverage_gr, gene_ranges)

  # Map to gene IDs
  coverage_df$gene_id <- NA_character_
  coverage_df$gene_id[queryHits(ov)] <- mcols(gene_ranges)$gene_id[subjectHits(ov)]

  if (!return_entrez) {
    # Convert Entrez IDs to gene SYMBOLs
    suppressMessages({coverage_df$gene_symbol <- mapIds(
      org.Hs.eg.db,
      keys      = gsub("^GeneID:", "", coverage_df$gene_id), # remove "GeneID:" prefix if present
      column    = "SYMBOL",
      keytype   = "ENTREZID",
      multiVals = "first"
    )})
  } else {
    coverage_df$gene_symbol <- coverage_df$gene_id
  }

  coverage_df
}


#' @export
# Function to compute cumulative genomic positions for coverage data and variant calls.
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
    summarize(chr_length = max(end)) %>%  # using max(end) as a proxy for length
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


#' @export
generateCNAPlotFromRegions2 <- function(depth_bigwig_file, variant_file, threshold = 1.5,
                                        downsample = 0.01, point_size = 0.5, colors = NULL,
                                        log_scale = TRUE, remove_centromeres = TRUE,
                                        max_value = NULL, min_value = NULL,
                                        min_variant_distance = 10000, samplename = "",
                                        show_cytobands = FALSE, cytoband_angle = 90,
                                        chr_filter = NULL, return_data = FALSE) {
  # Import coverage data from BigWig file and convert to data frame.
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- as.data.frame(coverage_data)

  # Filter for valid chromosomes (chr1 - chr22, chrX, chrY).
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))
  coverage_data <- coverage_data %>% filter(seqnames %in% valid_chromosomes)

  # Optionally, filter by a specific chromosome.
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


  # Read the VCF file and adjust column names.
  vcf <- fread(variant_file, skip = "#CHROM")
  colnames(vcf)[1] <- "CHROM"
  vcf <- vcf %>% filter(CHROM %in% valid_chromosomes)

  # Optionally, filter variant calls by chromosome.
  if (!is.null(chr_filter)) {
    vcf <- vcf %>% filter(CHROM == chr_filter)
  }

  # Extract variant information (end position and type from INFO field).
  end <- strsplit(vcf$INFO, ";") %>% sapply("[[", 3) %>% gsub("END=", "", .) %>% as.numeric()
  type <- strsplit(vcf$INFO, ";") %>% sapply("[[", 2) %>% gsub("SVTYPE=", "", .)
  variant_calls <- data.frame(seqnames = vcf$CHROM, start = vcf$POS, end = end, type = type)

  # Convert variant_calls to a GRanges object.
  variant_gr <- GRanges(
    seqnames = variant_calls$seqnames,
    ranges = IRanges(start = variant_calls$start, end = variant_calls$end)
  )

  # If cytoband annotation is desired, find overlaps with the cytogenetic bands.
  if (show_cytobands) {
    overlaps <- findOverlaps(variant_gr, hg19IdeogramCyto, type = "any")
    query_indices <- queryHits(overlaps)
    subject_indices <- subjectHits(overlaps)
    overlap_df <- data.frame(
      query = query_indices,
      cytoband = hg19IdeogramCyto$name[subject_indices]
    )

    # Use a standard if/else block inside summarize to avoid method-dispatch issues.
    cytoband_annotation <- overlap_df %>%
      mutate(cytoband = as.character(cytoband)) %>%  # ensure character type
      group_by(query) %>%
      summarize(
        cytoband_name = {
          if (n() == 1) {
            dplyr::first(cytoband)
          } else {
            paste(dplyr::first(cytoband), dplyr::last(cytoband), sep = "-")
          }
        }
      )

    variant_calls$cytoband <- cytoband_annotation$cytoband_name[match(1:nrow(variant_calls), cytoband_annotation$query)]
  } else {
    variant_calls$cytoband <- NA  # If not using cytoband labels, fill with NA.
  }

  # Filter variant calls by a minimum distance, if specified.
  if (!is.null(min_variant_distance) && min_variant_distance > 0) {
    variant_calls <- variant_calls %>%
      mutate(distance = end - start) %>%
      filter(distance >= min_variant_distance)
  }

  # Generate cumulative offsets for both coverage data and variant calls.
  offset <- generate_offsets(coverage_data = coverage_data, variant_calls = variant_calls)

  # Set default colors if not provided (the vector colors should be named appropriately).
  if (is.null(colors)) {
    chromosome_order <- c(as.character(1:22), "X", "Y")
    colors <- rep(c("black", "gray"), length.out = length(chromosome_order))
  }

  # Create a data frame for drawing variant highlight rectangles.
  if (nrow(offset$variant_calls) > 0) {
    data_highlight <- data.frame(
      xmin = offset$variant_calls$cumulative_start,
      xmax = offset$variant_calls$cumulative_end,
      ymin = rep(-Inf, nrow(offset$variant_calls)),
      ymax = rep(Inf, nrow(offset$variant_calls)),
      VariantType = offset$variant_calls$type,
      cytoband = offset$variant_calls$cytoband,
      seqnames = offset$variant_calls$seqnames,   # Chromosome names.
      x_center = (offset$variant_calls$cumulative_start + offset$variant_calls$cumulative_end) / 2
    )
    # Build a label that includes the chromosome and cytoband (e.g., "chr15: p11.1").
    data_highlight$label <- paste(data_highlight$seqnames, data_highlight$cytoband, sep = ": ")

    # Set text color based on the variant type.
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

  # Begin building the plot.
  p <- ggplot() +
    # Plot the coverage data points.
    geom_point(data = offset$coverage_data,
               aes(x = cumulative_genomic_position, y = delta, color = seqnames),
               size = point_size) +
    scale_color_manual(values = colors) +
    # Add variant highlight rectangles.
    geom_rect(data = data_highlight,
              aes(xmin = xmin, xmax = xmax,
                  ymin = ymin, ymax = ymax,
                  fill = VariantType),
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
    guides(color = "none")+
    scale_y_continuous(limits = c(ifelse(is.null(min_value), min(offset$coverage_data$delta, na.rm = TRUE), min_value),
                                  ifelse(is.null(max_value), max(offset$coverage_data$delta, na.rm = TRUE), max_value)))


  # Optionally add cytoband labels.
  if (show_cytobands) {
    # Choose a y position for labels that is within the defined y-axis range.
    label_y <- if (!is.null(max_value)) {
      max_value * 0.9
    } else {
      max(coverage_data$delta, na.rm = TRUE) * 0.9
    }

    # Exclude rows with missing cytoband values.
    data_highlight_labels <- data_highlight %>% filter(!is.na(cytoband))

    # Add the text labels with the specified rotation (cytoband_angle) and text colors.
    p <- p + geom_text(
      data = data_highlight_labels,
      aes(x = x_center, y = label_y, label = label, color = I(text_color)),
      size = 3, angle = cytoband_angle, hjust = 0
    )
  }

  # Optionally, return both the plot and the underlying data.
  if (return_data) {
    return(list(plot = p,
                coverage_data = offset$coverage_data,
                variant_calls = offset$variant_calls))
  } else {
    return(p)
  }
}
