library(ggplot2)
library(ggrepel)
library(GenomicRanges)
library(dplyr)
library(data.table)
library(rtracklayer)
library(org.Hs.eg.db)    # for gene symbol mapping, if needed

# We'll reuse your generate_offsets() from above

#' @export
generateCNAPlotDiscoverGenes2 <- function(
    depth_bigwig_file,
    variant_file,
    txdb,
    gene_delta_threshold = 2,
    downsample = 0.01,
    point_size = 0.5,
    colors = NULL,
    max_value = NULL,
    min_value = NULL,
    min_variant_distance = 10000,
    samplename = "",
    chr_filter = NULL,
    exclude_xy = T,   # New parameter to exclude chrX and chrY
    return_data = FALSE
) {
  # 1) Import coverage data
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- as.data.frame(coverage_data)

  # 2) Filter chromosomes
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))

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
  coverage_data$delta <- coverage_data$score - mean_coverage

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
      aes(x = cumulative_genomic_position, y = delta, color = seqnames),
      size = point_size
    ) +
    scale_color_manual(values = colors, drop = FALSE) +
    {
      if (nrow(data_highlight) > 0) {
        geom_rect(
          data = data_highlight,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = VariantType),
          color = "transparent",
          alpha = 0.2
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

#' @export
generateCNAPlotDiscoverGenes <- function(
    depth_bigwig_file,
    variant_file,
    txdb,
    gene_delta_threshold = 2,
    downsample = 0.01,
    point_size = 0.5,
    colors = NULL,
    max_value = NULL,
    min_value = NULL,
    min_variant_distance = 10000,
    samplename = "",
    chr_filter = NULL,
    exclude_xy = FALSE,
    highlight_genes = NULL,  # Vector of specific genes to label
    return_data = FALSE
) {
  # 1) Import coverage data
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- as.data.frame(coverage_data)

  # 2) Filter chromosomes
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))
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
  coverage_data$delta <- coverage_data$score - mean_coverage

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
      x_center    = (variant_df$cumulative_start + variant_df$cumulative_end) / 2
    )
  }

  # 9) Annotate genes (focus on 5' prime position)
  coverage_annot <- offset$coverage_data %>%
    mutate(end = start)  # Treat coverage as a single coordinate

  coverage_annot <- annotateCoverageWithGenes(coverage_annot, txdb = txdb)

  # Filter to highlight only specific genes
  if (!is.null(highlight_genes)) {
    coverage_hits <- coverage_annot %>%
      filter(gene_symbol %in% highlight_genes) %>%
      filter(abs(delta) >= gene_delta_threshold)
  } else {
    coverage_hits <- coverage_annot %>%
      filter(!is.na(gene_symbol)) %>%
      filter(abs(delta) >= gene_delta_threshold)
  }

  # Select the most 5' prime position for each gene
  coverage_hits <- coverage_hits %>%
    group_by(gene_symbol) %>%
    slice_min(start, n = 1) %>%  # Pick row with smallest start coordinate
    ungroup()

  # 10) Build base plot
  chromosome_midpoints <- offset$coverage_data %>%
    group_by(seqnames) %>%
    summarize(midpoint = mean(range(cumulative_genomic_position)), .groups = 'drop')

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
      aes(x = cumulative_genomic_position, y = delta, color = seqnames),
      size = point_size
    ) +
    scale_color_manual(values = colors, drop = FALSE) +
    {
      if (nrow(data_highlight) > 0) {
        geom_rect(
          data = data_highlight,
          aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = VariantType),
          color = "transparent",
          alpha = 0.2
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
      data = chromosome_midpoints,
      aes(xintercept = midpoint),
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
    scale_x_continuous(
      breaks = chromosome_midpoints$midpoint,
      labels = chromosome_midpoints$seqnames
    )

  # 11) Add gene labels at the most 5' prime position
  if (nrow(coverage_hits) > 0) {
    p <- p +
      geom_point(
        data = coverage_hits,
        aes(x = cumulative_genomic_position, y = delta),
        shape = 21,
        fill  = "white",
        color = "black",
        stroke = 1.2,
        size   = 3
      ) +
      geom_label_repel(
        data = coverage_hits,
        aes(
          x = cumulative_genomic_position,
          y = delta,
          label = gene_symbol
        ),
        box.padding      = 0.3,
        point.padding    = 0.3,
        min.segment.length = 0,
        nudge_y          = 50,
        color            = "black",
        fill             = "white",
        segment.color    = "black"
      )
  }

  # Return results
  if (return_data) {
    return(list(
      plot          = p,
      coverage_data = offset$coverage_data,
      variant_calls = offset$variant_calls,
      coverage_annot = coverage_annot
    ))
  } else {
    return(p)
  }
}
