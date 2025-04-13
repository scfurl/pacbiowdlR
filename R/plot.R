library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(dplyr)
library(data.table)
library(rtracklayer)
library(org.Hs.eg.db)    # for gene symbol mapping, if needed

# We'll reuse your generate_offsets() from above

#' @export
CNAPlot <- function(
    depth_bigwig_file,
    variant_file,
    txdb,
    method=c("fit", "delta", "loess"),
    gene_delta_threshold = 2,
    downsample = 0.01,
    point_size = 0.01,
    line_size = 0.1,
    colors = NULL,
    max_value = NULL,
    min_value = NULL,
    min_variant_distance = 10000,
    samplename = "",
    chr_filter = NULL,
    trend_window = 50,  # Number of points for running median trend line
    apply_weight = TRUE,  # If TRUE, weight delta values outside CNA calls
    outside_weight = 0.25,  # Multiplier for delta outside CNA calls
    inside_weight = 1,
    variant_alpha = 0.1,
    exclude_xy = T,   # New parameter to exclude chrX and chrY
    return_data = FALSE
) {
  method = method[1]


  # 1) Import coverage data
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- keepStandardChromosomes(coverage_data, pruning.mode = "tidy")
  coverage_data <- as.data.frame(coverage_data)

  # 2) Filter chromosomes
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))

  if(method=="fit"){
    data("gc_repeat_data")
    coverage_data$predicted_score <- gc_repeat_data$predicted_score
    coverage_data$repeat_fraction <- gc_repeat_data$repeat_fraction
    coverage_data$gc_content <- gc_repeat_data$gc_content
    coverage_data$width <- gc_repeat_data$width
  }

  # Exclude X and Y if requested
  if (exclude_xy) {
    valid_chromosomes <- paste0("chr", as.character(1:22))  # Remove X and Y
  }

  coverage_data <- coverage_data %>%
    filter(seqnames %in% valid_chromosomes)

  if (!is.null(chr_filter)) {
    coverage_data <- coverage_data %>% filter(seqnames == chr_filter)
  }

  # 3) Downsample coverage
  set.seed(42)
  sampled_indices <- sample(
    seq_len(nrow(coverage_data)),
    size = max(1, floor(nrow(coverage_data) * downsample))
  )
  coverage_data <- coverage_data[sampled_indices, ]

  # 4) Compute delta
  mean_coverage <- mean(coverage_data$score)
  if(method=="delta"){
    coverage_data$delta <- coverage_data$score - mean_coverage
  }

  if(method=="loess"){
    lcd <- coverage_data
    lo <- loess(score ~ start, data = lcd)
    lcd$loess<- lo$fitted
    lcd$repeat_fraction <- 1
    lcd <- compute_weighted_log2_ratios(lcd, observed_col = "score", reference_col = "loess", bin_length_col = "width")
    coverage_data$delta <- lcd$log2_ratio
    coverage_data$delta[is.infinite(coverage_data$delta)] <- 0
  }

  if(method=="fit"){
    mod <- lm(score~gc_content+repeat_fraction, data = coverage_data)
    coverage_data$predicted_score <- predict(mod, newdata = coverage_data)
    coverage_data <- compute_weighted_log2_ratios(coverage_data, observed_col = "score",
                                          reference_col = "predicted_score", bin_length_col = "width",
                                          repeat_fraction_col = "repeat_fraction")
    coverage_data$delta <- coverage_data$log2_ratio
    coverage_data$delta[is.infinite(coverage_data$delta)] <- 0
    if(!is.null(min_value)){
      coverage_data$delta[coverage_data$delta < min_value] <- min_value
    }
    if(!is.null(max_value)){
      coverage_data$delta[coverage_data$delta > max_value] <- max_value
    }
  }


  # 5) Read variants
  vcf <- fread(variant_file, skip = "#CHROM")
  colnames(vcf)[1] <- "CHROM"
  vcf <- vcf %>% filter(CHROM %in% valid_chromosomes)

  if (!is.null(chr_filter)) {
    vcf <- vcf %>% filter(CHROM == chr_filter)
  }

  end_vals <- strsplit(vcf$INFO, ";") %>%
    sapply("[[", 3) %>%
    gsub("END=", "", .) %>%
    as.numeric()

  type_vals <- strsplit(vcf$INFO, ";") %>%
    sapply("[[", 2) %>%
    gsub("SVTYPE=", "", .)

  variant_calls <- data.frame(
    seqnames = vcf$CHROM,
    start    = vcf$POS,
    end      = end_vals,
    type     = type_vals
  )

  # 6) Filter short variants
  if (!is.null(min_variant_distance) && min_variant_distance > 0) {
    variant_calls <- variant_calls %>%
      mutate(distance = end - start) %>%
      filter(distance >= min_variant_distance)
  }

  # 7) Compute offsets
  offset <- generate_offsets(coverage_data, variant_calls)

  # 8) Mark variants for rectangles
  variant_df <- offset$variant_calls
  data_highlight <- data.frame()
  if (nrow(variant_df) > 0) {
    data_highlight <- data.frame(
      xmin        = variant_df$cumulative_start,
      xmax        = variant_df$cumulative_end,
      ymin        = -Inf,
      ymax        = Inf,
      VariantType = variant_df$type,
      seqnames    = variant_df$seqnames,
      x_center    = (variant_df$cumulative_start + variant_df$cumulative_end)/2
    )
  }

  # 8b) OPTIONAL: Apply weighting to delta values
  # Create a weight multiplier for each point: if the point falls inside a CNA call, multiplier = 1;
  # otherwise, multiplier = outside_weight.
  offset$coverage_data <- offset$coverage_data %>%
    arrange(cumulative_genomic_position)

  if (apply_weight) {
    offset$coverage_data$multiplier <- sapply(offset$coverage_data$cumulative_genomic_position, function(pos) {
      if (nrow(data_highlight) > 0 && any(pos >= data_highlight$xmin & pos <= data_highlight$xmax))
        1
      else
        outside_weight
    })
  } else {
    offset$coverage_data$multiplier <- inside_weight
  }

  # Compute the weighted delta.
  offset$coverage_data$delta_weighted <- offset$coverage_data$delta * offset$coverage_data$multiplier

  # **************************
  # NEW: Compute running trend line over delta_weighted.
  # Use a rolling median (or change FUN to mean for a running mean).
  offset$coverage_data$trend <- rollapply(
    offset$coverage_data$delta_weighted,
    width = trend_window,
    FUN = median,
    fill = NA,
    align = "center"
  )
  # **************************


  # 9) Build base plot
  chromosome_boundaries <- offset$coverage_data %>%
    group_by(seqnames) %>%
    summarize(boundary = min(cumulative_genomic_position), .groups = 'drop')

  chromosome_midpoints <- offset$coverage_data %>%
    group_by(seqnames) %>%
    summarize(midpoint = mean(range(cumulative_genomic_position)), .groups = 'drop')

  # If user didn't supply color vector, just alternate black/gray
  if (is.null(colors)) {
    chromosome_order <- c(as.character(1:22), "X", "Y")

    if (exclude_xy) {
      chromosome_order <- as.character(1:22)  # Remove X and Y
    }

    colors <- rep(c("black", "gray"), length.out = length(chromosome_order))
    names(colors) <- paste0("chr", chromosome_order)
  }

  p <- ggplot() +
    geom_point(
      data = offset$coverage_data,
      aes(x = cumulative_genomic_position, y = delta_weighted, color = seqnames),
      size = point_size
    ) +
    scale_color_manual(values = colors, drop = FALSE) +
    # Add the trend line from the running median
    geom_line(
      data = offset$coverage_data,
      aes(x = cumulative_genomic_position, y = trend),
      color = "black", size = line_size, na.rm = TRUE
    ) +
    {
      if (nrow(data_highlight) > 0) {
        geom_rect(
          data = data_highlight,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = VariantType),
          color = "transparent",
          alpha = variant_alpha
        )
      } else {
        NULL
      }
    } +
    scale_fill_manual(
      values = c("DEL" = "blue", "DUP" = "red"),
      na.value = "grey80"
    ) +
    geom_vline(
      data = chromosome_boundaries,
      aes(xintercept = boundary),
      color = "grey",
      linetype = "dashed",
      alpha = 0.6
    ) +
    labs(
      title = paste0("Genome-Wide CNA ", samplename),
      x = "Chromosome",
      y = "Delta"
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.ticks.x = element_blank(),
      panel.grid   = element_blank()
    ) +
    guides(color = "none") +
    scale_y_continuous(
      limits = c(
        ifelse(is.null(min_value), min(offset$coverage_data$delta, na.rm = TRUE), min_value),
        ifelse(is.null(max_value), max(offset$coverage_data$delta, na.rm = TRUE), max_value)
      )
    ) +
    scale_x_continuous(
      breaks = chromosome_midpoints$midpoint,
      labels = chromosome_midpoints$seqnames
    )

  # Return results
  if (return_data) {
    return(list(
      plot          = p,
      coverage_data = offset$coverage_data,
      variant_calls = offset$variant_calls
    ))
  } else {
    return(p)
  }
}
library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(dplyr)
library(data.table)
library(rtracklayer)
library(org.Hs.eg.db)    # for gene symbol mapping, if needed
library(zoo)             # for rollapply

#' @export
CNAPlot_Highlight <- function(
    depth_bigwig_file,
    variant_file,
    txdb,
    gene_delta_threshold = 2,
    downsample = 0.01,
    point_size = 0.01,
    line_size = 0.1,
    colors = NULL,
    max_value = NULL,
    min_value = NULL,
    min_variant_distance = 10000,
    method = c("fit", "delta", "loess"),
    trend_window = 50,         # Number of points for running median trend line
    apply_weight = TRUE,       # If TRUE, weight delta values outside CNA calls
    outside_weight = 0.25,     # Multiplier for delta outside CNA calls
    inside_weight = 1,
    variant_alpha = 0.1,
    nudge_y = 3,
    samplename = "",
    chr_filter = NULL,
    exclude_xy = FALSE,
    highlight_genes = NULL,    # Vector of specific genes to label
    gband_file = file.path("~/develop/pacbiowdlR/cytobands.tsv"),         # NEW: file path to cytogenetic band data
    gband_y_offset = -0.2,      # Vertical offset for G-band labels (in y units)
    gband_text_size = 2,       # Text size for G-band labels
    showCNAbands = F,
    return_data = FALSE
) {
  method <- method[1]

  # 1) Import coverage data
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- keepStandardChromosomes(coverage_data, pruning.mode = "tidy")
  coverage_data <- as.data.frame(coverage_data)

  # 2) Filter chromosomes
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))
  if(method == "fit"){
    data("gc_repeat_data")
    coverage_data$predicted_score <- gc_repeat_data$predicted_score
    coverage_data$repeat_fraction <- gc_repeat_data$repeat_fraction
    coverage_data$gc_content <- gc_repeat_data$gc_content
    coverage_data$width <- gc_repeat_data$width
  }
  if (exclude_xy) {
    valid_chromosomes <- paste0("chr", as.character(1:22))
  }
  coverage_data <- coverage_data %>% filter(seqnames %in% valid_chromosomes)
  if (!is.null(chr_filter)) {
    coverage_data <- coverage_data %>% filter(seqnames == chr_filter)
  }

  # 3) Downsample coverage data
  set.seed(42)
  sampled_indices <- sample(seq_len(nrow(coverage_data)),
                            size = max(1, floor(nrow(coverage_data) * downsample)))
  coverage_data <- coverage_data[sampled_indices, ]

  # 4) Compute delta values
  mean_coverage <- mean(coverage_data$score)
  if(method == "delta"){
    coverage_data$delta <- coverage_data$score - mean_coverage
  } else if(method == "loess"){
    lcd <- coverage_data
    lo <- loess(score ~ start, data = lcd)
    lcd$loess <- lo$fitted
    lcd$repeat_fraction <- 1
    lcd <- compute_weighted_log2_ratios(lcd,
                                        observed_col = "score",
                                        reference_col = "loess",
                                        bin_length_col = "width")
    coverage_data$delta <- lcd$log2_ratio
    coverage_data$delta[is.infinite(coverage_data$delta)] <- 0
  } else if(method == "fit"){
    mod <- lm(score ~ gc_content + repeat_fraction, data = coverage_data)
    coverage_data$predicted_score <- predict(mod, newdata = coverage_data)
    coverage_data <- compute_weighted_log2_ratios(coverage_data,
                                                  observed_col = "score",
                                                  reference_col = "predicted_score",
                                                  bin_length_col = "width",
                                                  repeat_fraction_col = "repeat_fraction")
    coverage_data$delta <- coverage_data$log2_ratio
    coverage_data$delta[is.infinite(coverage_data$delta)] <- 0
    if(!is.null(min_value)){
      coverage_data$delta[coverage_data$delta < min_value] <- min_value
    }
    if(!is.null(max_value)){
      coverage_data$delta[coverage_data$delta > max_value] <- max_value
    }
  }

  # 5) Read variants from VCF
  vcf <- fread(variant_file, skip = "#CHROM")
  colnames(vcf)[1] <- "CHROM"
  vcf <- vcf %>% filter(CHROM %in% valid_chromosomes)
  if (!is.null(chr_filter)) {
    vcf <- vcf %>% filter(CHROM == chr_filter)
  }
  end_vals <- strsplit(vcf$INFO, ";") %>% sapply("[[", 3) %>% gsub("END=", "", .) %>% as.numeric()
  type_vals <- strsplit(vcf$INFO, ";") %>% sapply("[[", 2) %>% gsub("SVTYPE=", "", .)
  variant_calls <- data.frame(
    seqnames = vcf$CHROM,
    start    = vcf$POS,
    end      = end_vals,
    type     = type_vals
  )

  # 6) Filter short variants
  if (!is.null(min_variant_distance) && min_variant_distance > 0) {
    variant_calls <- variant_calls %>% mutate(distance = end - start) %>% filter(distance >= min_variant_distance)
  }

  # 7) Compute offsets (assumes generate_offsets() is defined elsewhere)
  offset <- generate_offsets(coverage_data, variant_calls)

  # 8) Mark variant regions for rectangles
  variant_df <- offset$variant_calls
  data_highlight <- data.frame()
  if (nrow(variant_df) > 0) {
    data_highlight <- data.frame(
      xmin     = variant_df$cumulative_start,
      xmax     = variant_df$cumulative_end,
      ymin     = -Inf,
      ymax     = Inf,
      VariantType = variant_df$type,
      seqnames = variant_df$seqnames,
      x_center = (variant_df$cumulative_start + variant_df$cumulative_end) / 2
    )
  }

  # 8b) Apply weighting to delta values
  offset$coverage_data <- offset$coverage_data %>% arrange(cumulative_genomic_position)
  if (apply_weight) {
    offset$coverage_data$multiplier <- sapply(offset$coverage_data$cumulative_genomic_position, function(pos) {
      if (nrow(data_highlight) > 0 && any(pos >= data_highlight$xmin & pos <= data_highlight$xmax))
        1
      else
        outside_weight
    })
  } else {
    offset$coverage_data$multiplier <- inside_weight
  }
  offset$coverage_data$delta_weighted <- offset$coverage_data$delta * offset$coverage_data$multiplier

  # 8c) Compute running trend line over weighted delta values
  offset$coverage_data$trend <- rollapply(
    offset$coverage_data$delta_weighted,
    width = trend_window,
    FUN = median,   # use mean here if you prefer a running mean instead
    fill = NA,
    align = "center"
  )

  # 9) Annotate genes (focus on 5' prime position)
  coverage_annot <- offset$coverage_data %>% mutate(end = start)  # treat as a single coordinate
  # Use single.strand.genes.only = FALSE so that multi-strand genes aren’t dropped
  coverage_annot <- annotateCoverageWithGenes(coverage_annot, txdb = txdb)

  # For gene labels, if highlight_genes is provided, retain those genes regardless of delta threshold;
  # otherwise, filter by gene_delta_threshold.
  if (!is.null(highlight_genes)) {
    coverage_hits <- coverage_annot %>% filter(gene_symbol %in% highlight_genes)
  } else {
    coverage_hits <- coverage_annot %>% filter(!is.na(gene_symbol)) %>% filter(abs(delta) >= gene_delta_threshold)
  }

  # Select one representative coordinate (the most 5'-prime) per gene
  coverage_hits <- coverage_hits %>% group_by(gene_symbol) %>% slice_min(start, n = 1) %>% ungroup()

  # 10) Build base plot (get chromosome midpoints for x-axis labels)
  chromosome_midpoints <- offset$coverage_data %>% group_by(seqnames) %>% summarize(midpoint = mean(range(cumulative_genomic_position)), .groups = 'drop')
  if (is.null(colors)) {
    chromosome_order <- c(as.character(1:22), "X", "Y")
    if (exclude_xy) {
      chromosome_order <- as.character(1:22)
    }
    colors <- rep(c("black", "gray"), length.out = length(chromosome_order))
    names(colors) <- paste0("chr", chromosome_order)
  }

  p <- ggplot() +
    geom_point(
      data = offset$coverage_data,
      aes(x = cumulative_genomic_position, y = delta_weighted, color = seqnames),
      size = point_size
    ) +
    geom_line(
      data = offset$coverage_data,
      aes(x = cumulative_genomic_position, y = trend),
      color = "black", size = line_size, na.rm = TRUE
    ) +
    scale_color_manual(values = colors, drop = FALSE) +
    {
      if (nrow(data_highlight) > 0) {
        geom_rect(
          data = data_highlight,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = VariantType),
          color = "transparent", alpha = variant_alpha
        )
      } else {
        NULL
      }
    } +
    scale_fill_manual(
      values = c("DEL" = "blue", "DUP" = "red"),
      na.value = "grey80"
    ) +
    labs(
      title = paste0("Genome-Wide CNA ", samplename),
      x = "Chromosome",
      y = "Delta"
    ) +
    theme_minimal() +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1),
      axis.ticks.x = element_blank(),
      panel.grid   = element_blank()
    ) +
    guides(color = "none") +
    scale_x_continuous(
      breaks = chromosome_midpoints$midpoint,
      labels = chromosome_midpoints$seqnames
    ) +
    scale_y_continuous(
      limits = c(
        ifelse(is.null(min_value), min(offset$coverage_data$delta_weighted, na.rm = TRUE), min_value),
        ifelse(is.null(max_value), max(offset$coverage_data$delta_weighted, na.rm = TRUE), max_value)
      )
    )

  # 11) Add gene labels at the most 5' prime position
  if (nrow(coverage_hits) > 0) {
    p <- p +
      geom_point(
        data = coverage_hits,
        aes(x = cumulative_genomic_position, y = delta_weighted),
        shape = 21, fill = "white", color = "black", stroke = 1.2, size = 3
      ) +
      geom_label_repel(
        data = coverage_hits,
        aes(x = cumulative_genomic_position, y = delta_weighted, label = gene_symbol),
        box.padding = 0.3, point.padding = 0.3, min.segment.length = 0,
        nudge_y = nudge_y, color = "black", fill = "white", segment.color = "black"
      )
  }

  # 12) NEW: Annotate G-band identity for variant regions (from the VCF)
  # This step uses an external band file if provided.
  if (!is.null(gband_file) & showCNAbands) {
    # Read the band file (assumed to be tab-delimited with columns: contig, start, end, name, giemsa)
    bands_df <- fread(gband_file)
    # Create a 'chr' column (e.g., "chr1") from contig
    bands_df$chr <- paste0("chr", bands_df$contig)
    band_gr <- GRanges(seqnames = bands_df$chr,
                       ranges = IRanges(start = bands_df$start, end = bands_df$end),
                       band = bands_df$name,
                       giemsa = bands_df$giemsa)
    bands_gr <- keepStandardChromosomes(band_gr, pruning.mode = "tidy")
    # Convert variant_calls into GRanges (assumed to be in the same order as data_highlight)
    var_gr <- GRanges(seqnames = variant_calls$seqnames,
                      ranges = IRanges(start = variant_calls$start, end = variant_calls$end))

    # Find overlaps between variant calls and bands
    ov <- findOverlaps(var_gr, band_gr)
    if (length(ov) > 0) {
      # For each variant, take the first overlapping band (or you can refine this selection)
      var_band <- data.frame(
        variant_idx = queryHits(ov),
        band = mcols(band_gr)$band[subjectHits(ov)],
        giemsa = mcols(band_gr)$giemsa[subjectHits(ov)]
      )
      var_band <- var_band %>% group_by(variant_idx) %>% slice(1) %>% ungroup()

      # Add band annotation to the variant_calls and then to data_highlight
      variant_calls$band <- NA
      variant_calls$giemsa <- NA
      variant_calls$band[var_band$variant_idx] <- var_band$band
      variant_calls$giemsa[var_band$variant_idx] <- var_band$giemsa

      data_highlight$band <- variant_calls$band
      data_highlight$giemsa <- variant_calls$giemsa

      # For band annotation, create a data frame with label text and x position
      band_annotations <- data_highlight %>%
        filter(!is.na(band)) %>%
        mutate(label = band,
               x_label = (xmin + xmax)/2)

      # Compute y position for labels (place these just at the top)
      y_max <- ifelse(is.null(max_value), min(offset$coverage_data$delta_weighted, na.rm = TRUE), max_value)
      y_label <- y_max + gband_y_offset

      p <- p +
        geom_text(
          data = band_annotations,
          aes(x = x_label, label = label, color = VariantType),
          y = y_label,
          size = gband_text_size,
          vjust = 1
        )
    }
  }

  # 13) Return results
  if (return_data) {
    return(list(
      plot = p,
      coverage_data = offset$coverage_data,
      variant_calls = offset$variant_calls,
      coverage_annot = coverage_annot
    ))
  } else {
    return(p)
  }
}

# Helper functions ---------------------------------------------------------

compute_weighted_log2_ratios <- function(data,
                                         observed_col = "observed",
                                         reference_col = "reference",
                                         bin_length_col = "bin_length",
                                         repeat_fraction_col = "repeat_fraction") {
  required_cols <- c(observed_col, reference_col, bin_length_col, repeat_fraction_col)
  if (!all(required_cols %in% colnames(data))) {
    stop("One or more required columns are missing from the data.")
  }
  data$log2_ratio <- with(data, log2(get(observed_col) / get(reference_col)))
  median_length <- median(data[[bin_length_col]], na.rm = TRUE)
  data$weight <- with(data,
                      (get(bin_length_col) / median_length) * (1 - get(repeat_fraction_col)))
  return(data)
}

weighted_segment_log2 <- function(segment_data) {
  total_weight <- sum(segment_data$weight, na.rm = TRUE)
  if (total_weight == 0) return(NA)
  weighted_avg <- sum(segment_data$log2_ratio * segment_data$weight, na.rm = TRUE) / total_weight
  return(weighted_avg)
}
