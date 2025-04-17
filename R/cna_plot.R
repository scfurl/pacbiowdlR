#' CNA Plot: Generate a Genome-Wide CNA Plot
#'
#' This function creates a genome-wide copy number alteration (CNA) plot from a BigWig depth file and a VCF variant file.
#' It computes delta values using one of three methods ("fit", "delta", or "loess"), downsampling the coverage data,
#' and overlays variant calls as rectangles. A running median trend line is computed on weighted delta values.
#'
#' @param depth_bigwig_file Character. File path to a BigWig file containing depth/coverage information.
#' @param variant_file Character. File path to a variant call file (VCF) containing structural variant information.
#' @param txdb A transcript database object (e.g., from the GenomicFeatures package) used for gene annotation.
#' @param method Character. One of \code{"fit"}, \code{"delta"}, or \code{"loess"} determining the method to compute delta values.
#'   Only the first element of the provided vector is used.
#' @param gene_delta_threshold Numeric. Threshold applied to the delta values for gene annotation.
#' @param downsample Numeric. Proportion of coverage data to retain (e.g., 0.01 for 1 percent).
#' @param point_size Numeric. Size of the plotted data points.
#' @param line_size Numeric. Size of the plotted trend line.
#' @param line_color Character. Color for the trend line.
#' @param colors Named vector of colors for chromosomes. If \code{NULL}, an alternating palette of black/gray is used.
#' @param max_value Numeric. Maximum allowed delta value; values above this are capped.
#' @param min_value Numeric. Minimum allowed delta value; values below this are capped.
#' @param min_variant_distance Numeric. Minimum distance (in bp) for a variant call to be retained.
#' @param samplename Character. A label for the sample that is appended to the plot title.
#' @param chr_filter Character. If specified, only data from the given chromosome are processed.
#' @param trend_window Integer. The number of consecutive data points over which to compute the running median trend line.
#' @param apply_weight Logical. If \code{TRUE}, weight multipliers are applied to delta values outside CNA calls.
#' @param outside_weight Numeric. Multiplier applied to delta values outside CNA calls.
#' @param inside_weight Numeric. Multiplier applied to delta values inside CNA calls.
#' @param variant_alpha Numeric. Transparency level for variant rectangles.
#' @param trend_regions Character. One of \code{"both"}, \code{"inside"}, or \code{"outside"} to control where the trend line is plotted.
#'   Default is \code{"both"}.
#' @param exclude_xy Logical. If \code{TRUE}, chromosomes X and Y are excluded from the analysis.
#' @param return_data Logical. If \code{TRUE}, returns a list containing the plot, the coverage data, and variant calls.
#'
#' @return Either a \code{ggplot2} object representing the CNA plot or a list with additional processed data when \code{return_data} is \code{TRUE}.
#'
#' @details The function performs several steps:
#' \enumerate{
#'   \item Imports coverage data from the BigWig file using \code{import.bw()} and keeps only standard chromosomes (with \code{keepStandardChromosomes()}).
#'   \item Optionally loads external GC/repeat data when using the "fit" method.
#'   \item Filters and down-samples the coverage data.
#'   \item Computes delta values using the specified method:
#'         \itemize{
#'           \item[\code{"delta"}] subtracts the mean coverage.
#'           \item[\code{"loess"}] fits a LOESS model and computes log2 ratios.
#'           \item[\code{"fit"}] fits a linear model to predict coverage based on GC content and repeat fraction.
#'         }
#'   \item Reads and filters variant calls from the VCF file.
#'   \item Computes genomic offsets via \code{generate_offsets()} (assumed to be defined elsewhere).
#'   \item Applies weighting to delta values and computes a running median trend line using \code{rollapply()}.
#'   \item Based on the value of \code{trend_regions}, the trend line is shown only for regions:
#'         \itemize{
#'           \item \code{"inside"}: only within CNV call regions,
#'           \item \code{"outside"}: only outside CNV call regions, or
#'           \item \code{"both"}: across all regions (the default).
#'         }
#'   \item Builds the plot with \code{ggplot2} incorporating points, trend line, variant rectangles, and chromosome boundaries.
#' }
#'
#' @note This function assumes that the helper function \code{generate_offsets()} is defined in your package.
#'
#' @importFrom rtracklayer import.bw
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom data.table fread
#' @importFrom dplyr filter arrange group_by summarize mutate slice_min ungroup
#' @importFrom zoo rollapply
#' @importFrom ggplot2 ggplot geom_point geom_line scale_color_manual scale_fill_manual geom_vline labs theme_minimal theme guides scale_y_continuous scale_x_continuous element_text element_blank aes
#' @importFrom ggrepel geom_label_repel
#' @importFrom magrittr %>%
#' @keywords internal
#' @export
CNA_plot <- function(
    depth_bigwig_file,
    variant_file,
    txdb,
    method = c("fit", "delta", "loess"),
    gene_delta_threshold = 2,
    downsample = 0.1,
    point_size = 0.01,
    line_size = 0.1,
    line_color = "red",
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
    trend_regions = "inside",
    exclude_xy = TRUE,
    return_data = FALSE
) {
  method <- method[1]

  # 1) Import coverage data using rtracklayer
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- keepStandardChromosomes(coverage_data, pruning.mode = "tidy")
  coverage_data <- as.data.frame(coverage_data)

  # 2) Filter chromosomes: standard chromosomes chr1-22, X, Y.
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))
  if (method == "fit") {
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
  sampled_indices <- sample(
    seq_len(nrow(coverage_data)),
    size = max(1, floor(nrow(coverage_data) * downsample))
  )
  coverage_data <- coverage_data[sampled_indices, ]

  # 4) Compute delta values
  mean_coverage <- mean(coverage_data$score)
  if (method == "delta") {
    coverage_data$delta <- coverage_data$score - mean_coverage
  }
  if (method == "loess") {
    lcd <- coverage_data
    lo <- loess(score ~ start, data = lcd)
    lcd$loess <- lo$fitted
    lcd$repeat_fraction <- 1
    lcd <- compute_weighted_log2_ratios(lcd, observed_col = "score",
                                        reference_col = "loess",
                                        bin_length_col = "width")
    coverage_data$delta <- lcd$log2_ratio
    coverage_data$delta[is.infinite(coverage_data$delta)] <- 0
  }
  if (method == "fit") {
    mod <- lm(score ~ gc_content + repeat_fraction, data = coverage_data)
    coverage_data$predicted_score <- predict(mod, newdata = coverage_data)
    coverage_data <- compute_weighted_log2_ratios(coverage_data, observed_col = "score",
                                                  reference_col = "predicted_score",
                                                  bin_length_col = "width",
                                                  repeat_fraction_col = "repeat_fraction")
    coverage_data$delta <- coverage_data$log2_ratio
    coverage_data$delta[is.infinite(coverage_data$delta)] <- 0
  }

  # 5) Read variants from file
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
    variant_calls <- variant_calls %>%
      mutate(distance = end - start) %>%
      filter(distance >= min_variant_distance)
  }

  # 7) Compute offsets (assumes generate_offsets() is defined elsewhere)
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

  # 8b) Optional: Apply weighting to delta values
  offset$coverage_data <- offset$coverage_data %>% arrange(cumulative_genomic_position)
  if (apply_weight) {
    offset$coverage_data$multiplier <- sapply(offset$coverage_data$cumulative_genomic_position, function(pos) {
      if (nrow(data_highlight) > 0 && any(pos >= data_highlight$xmin & pos <= data_highlight$xmax))
        inside_weight
      else
        outside_weight
    })
  } else {
    offset$coverage_data$multiplier <- inside_weight
  }
  offset$coverage_data$delta_weighted <- offset$coverage_data$delta * offset$coverage_data$multiplier

  # 8c) Compute running trend line (using a rolling median)
  offset$coverage_data$trend <- rollapply(
    offset$coverage_data$delta_weighted,
    width = trend_window,
    FUN = median,
    fill = NA,
    align = "center"
  )

  # --- NEW: Subset or mask the trend line by region ----
  # Create a logical flag vector that is TRUE when a point is inside any variant call region
  inside_flag <- sapply(offset$coverage_data$cumulative_genomic_position, function(pos) {
    if (nrow(data_highlight) > 0 && any(pos >= data_highlight$xmin & pos <= data_highlight$xmax)) {
      TRUE
    } else {
      FALSE
    }
  })
  # Adjust the trend line based on the 'trend_regions' parameter:
  # "both": no change, "inside": set trend to NA for points outside variant regions,
  # "outside": set trend to NA for points inside variant regions.
  if (trend_regions == "inside") {
    offset$coverage_data$trend[!inside_flag] <- NA
  } else if (trend_regions == "outside") {
    offset$coverage_data$trend[inside_flag] <- NA
  }

  # 9) Build base plot
  # Apply min/max constraints to the weighted values (if provided)
  if (!is.null(min_value)) {
    offset$coverage_data$delta_weighted[offset$coverage_data$delta_weighted < min_value] <- min_value
  }
  if (!is.null(max_value)) {
    offset$coverage_data$delta_weighted[offset$coverage_data$delta_weighted > max_value] <- max_value
  }
  chromosome_boundaries <- offset$coverage_data %>%
    group_by(seqnames) %>%
    summarize(boundary = min(cumulative_genomic_position), .groups = 'drop')
  chromosome_midpoints <- offset$coverage_data %>%
    group_by(seqnames) %>%
    summarize(midpoint = mean(range(cumulative_genomic_position)), .groups = 'drop')

  # If colors are not supplied, alternate black/gray
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
    scale_color_manual(values = colors, drop = FALSE) +
    # Add the running trend line (only the non-NA values will be plotted)
    geom_line(
      data = offset$coverage_data,
      aes(x = cumulative_genomic_position, y = trend),
      color = line_color, size = line_size, na.rm = TRUE
    ) +
    { if (nrow(data_highlight) > 0) {
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
      title = paste0("CNA ", samplename),
      x = "Chromosome",
      y = ifelse(method %in% c("loess", "fit"),
                 expression(log[2](frac(Obs, Exp)) ~ "Coverage"),
                 expression(delta ~ Coverage))
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
        ifelse(is.null(min_value), min(offset$coverage_data$delta_weighted, na.rm = TRUE), min_value),
        ifelse(is.null(max_value), max(offset$coverage_data$delta_weighted, na.rm = TRUE), max_value)
      )
    ) +
    scale_x_continuous(
      breaks = chromosome_midpoints$midpoint,
      labels = chromosome_midpoints$seqnames
    )

  # Return either the plot or a list with additional data
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

#' CNA Plot with Gene Highlighting
#'
#' Generates a genome-wide CNA plot with additional gene annotation and labeling.
#' This function extends \code{CNA_plot} by annotating gene coordinates (using
#' \code{annotateCoverageWithGenes()}) and optionally highlighting selected genes
#' and cytogenetic (G-band) features.
#'
#' @param depth_bigwig_file Character. File path to a BigWig file containing coverage depth data.
#' @param variant_file Character. File path to a VCF file with variant data.
#' @param txdb A transcript database object used for gene annotation.
#' @param org A species annotation object (e.g., \code{org.Hs.eg.db}) for mapping gene identifiers.
#' @param gene_delta_threshold Numeric. Delta threshold for gene annotation filtering.
#' @param downsample Numeric. Proportion of the coverage data to retain after downsampling.
#' @param point_size Numeric. Size of individual points in the plot.
#' @param line_size Numeric. Size of the trend line.
#' @param line_color Character. Color of the trend line.
#' @param colors Named vector of colors for chromosomes; if \code{NULL}, defaults to an alternating palette.
#' @param max_value Numeric. Maximum allowed delta value; values above this are capped.
#' @param min_value Numeric. Minimum allowed delta value; values below this are capped.
#' @param min_variant_distance Numeric. Minimum variant length (in bp) to include.
#' @param method Character. One of \code{"fit"}, \code{"delta"}, or \code{"loess"} determining the method used to compute delta values.
#' @param trend_window Integer. Window size for computing the running median trend line.
#' @param apply_weight Logical. If \code{TRUE}, applies a weight multiplier to delta values outside CNA calls.
#' @param outside_weight Numeric. Weight multiplier applied to delta values outside CNA calls.
#' @param inside_weight Numeric. Weight multiplier applied to delta values inside CNA calls.
#' @param variant_alpha Numeric. Transparency level for variant call rectangles.
#' @param nudge_y Numeric. Vertical nudge for positioning gene labels.
#' @param samplename Character. Sample label used in the plot title.
#' @param chr_filter Character. If provided, restricts the analysis to the specified chromosome.
#' @param exclude_xy Logical. If \code{TRUE}, chromosomes X and Y are excluded.
#' @param highlight_genes Character vector. Specific gene symbols to highlight in the plot.
#' @param gband_file Character. File path to a cytogenetic band data file (e.g., a TSV file).
#' @param gband_y_offset Numeric. Vertical offset for G-band label placement.
#' @param gband_text_size Numeric. Text size for G-band labels.
#' @param showCNAbands Logical. If \code{TRUE}, displays G-band annotations on the plot.
#' @param trend_regions Character. One of \code{"both"}, \code{"inside"}, or \code{"outside"}.
#'   This option controls where the running median trend line is displayed:
#'   \describe{
#'     \item{\code{"both"}}{Trend line is plotted for all points (default).}
#'     \item{\code{"inside"}}{Trend line is only plotted for positions that fall within variant regions.}
#'     \item{\code{"outside"}}{Trend line is only plotted for positions that fall outside variant regions.}
#'   }
#' @param return_data Logical. If \code{TRUE}, returns a list containing the plot, the coverage data, variant calls, and gene annotations.
#'
#' @return Either a \code{ggplot2} object representing the CNA plot with gene highlights, or a list with additional data when \code{return_data} is \code{TRUE}.
#'
#' @details The function performs several steps:
#' \enumerate{
#'   \item Imports and filters coverage data from the BigWig file using \code{import.bw()} and retains standard chromosomes (using \code{keepStandardChromosomes()}).
#'   \item Optionally loads external GC/repeat data when using the \code{"fit"} method.
#'   \item Downsamples the coverage data.
#'   \item Computes delta values using the specified method:
#'         \itemize{
#'           \item[\code{"delta"}] subtracts the mean coverage.
#'           \item[\code{"loess"}] fits a LOESS model and computes log2 ratios.
#'           \item[\code{"fit"}] fits a linear model to predict coverage based on GC content and repeat fraction.
#'         }
#'   \item Reads and filters variant calls from the VCF file.
#'   \item Computes genomic offsets via \code{generate_offsets()} and marks variant regions.
#'   \item Applies weighting to delta values with different multipliers for points inside and outside CNV calls.
#'   \item Computes a running median trend line using \code{rollapply()}, and, based on \code{trend_regions},
#'         subsets the trend line to display only points inside variants, outside variants, or in both regions.
#'   \item Annotates gene information on the coverage data using \code{annotateCoverageWithGenes()}.
#'   \item Constructs the final plot with \code{ggplot2} showing data points, the trend line, variant rectangles, and gene labels (with optional G-band annotations).
#' }
#'
#' @importFrom rtracklayer import.bw
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom data.table fread
#' @importFrom dplyr filter arrange group_by summarize mutate slice_min ungroup
#' @importFrom zoo rollapply
#' @importFrom ggplot2 ggplot geom_point geom_line scale_color_manual scale_fill_manual geom_vline labs theme_minimal theme guides scale_x_continuous scale_y_continuous element_text element_blank geom_text
#' @importFrom ggrepel geom_label_repel
#' @importFrom S4Vectors subjectHits queryHits
#' @importFrom magrittr "%>%"
#'
#' @keywords internal
#' @export
CNA_plot_highlight <- function(
    depth_bigwig_file,
    variant_file,
    txdb,
    org,
    gene_delta_threshold = 2,
    downsample = 0.1,
    point_size = 0.01,
    line_size = 0.1,
    line_color = "red",
    colors = NULL,
    max_value = NULL,
    min_value = NULL,
    min_variant_distance = 10000,
    method = c("fit", "delta", "loess"),
    trend_window = 50,
    apply_weight = TRUE,
    outside_weight = 0.25,
    inside_weight = 1,
    variant_alpha = 0.1,
    nudge_y = 3,
    samplename = "",
    chr_filter = NULL,
    exclude_xy = FALSE,
    highlight_genes = NULL,
    gband_file = file.path("~/develop/pacbiowdlR/cytobands.tsv"),
    gband_y_offset = -0.2,
    gband_text_size = 2,
    showCNAbands = FALSE,
    trend_regions = "both",
    return_data = FALSE
) {
  method <- method[1]

  # 1) Import coverage data using rtracklayer
  coverage_data <- import.bw(depth_bigwig_file)
  coverage_data <- keepStandardChromosomes(coverage_data, pruning.mode = "tidy")
  coverage_data <- as.data.frame(coverage_data)

  # 2) Filter chromosomes: standard chromosomes chr1-22, X, Y.
  valid_chromosomes <- paste0("chr", c(as.character(1:22), "X", "Y"))
  if (method == "fit") {
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

  # 4) Compute delta values using the specified method
  mean_coverage <- mean(coverage_data$score)
  if (method == "delta") {
    coverage_data$delta <- coverage_data$score - mean_coverage
  } else if (method == "loess") {
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
  } else if (method == "fit") {
    mod <- lm(score ~ gc_content + repeat_fraction, data = coverage_data)
    coverage_data$predicted_score <- predict(mod, newdata = coverage_data)
    coverage_data <- compute_weighted_log2_ratios(coverage_data,
                                                  observed_col = "score",
                                                  reference_col = "predicted_score",
                                                  bin_length_col = "width",
                                                  repeat_fraction_col = "repeat_fraction")
    coverage_data$delta <- coverage_data$log2_ratio
    coverage_data$delta[is.infinite(coverage_data$delta)] <- 0
  }

  # 5) Read variants from file
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
      xmin        = variant_df$cumulative_start,
      xmax        = variant_df$cumulative_end,
      ymin        = -Inf,
      ymax        = Inf,
      VariantType = variant_df$type,
      seqnames    = variant_df$seqnames,
      x_center    = (variant_df$cumulative_start + variant_df$cumulative_end) / 2
    )
  }

  # 8b) Optional: Apply weighting to delta values
  offset$coverage_data <- offset$coverage_data %>% arrange(cumulative_genomic_position)
  if (apply_weight) {
    offset$coverage_data$multiplier <- sapply(offset$coverage_data$cumulative_genomic_position, function(pos) {
      if (nrow(data_highlight) > 0 && any(pos >= data_highlight$xmin & pos <= data_highlight$xmax))
        inside_weight
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
    FUN = median,
    fill = NA,
    align = "center"
  )

  # --- Subset the trend line according to the trend_regions parameter ---
  # Create a flag for points that fall within any variant region
  inside_flag <- sapply(offset$coverage_data$cumulative_genomic_position, function(pos) {
    if (nrow(data_highlight) > 0 && any(pos >= data_highlight$xmin & pos <= data_highlight$xmax)) {
      TRUE
    } else {
      FALSE
    }
  })

  if (trend_regions == "inside") {
    offset$coverage_data$trend[!inside_flag] <- NA
  } else if (trend_regions == "outside") {
    offset$coverage_data$trend[inside_flag] <- NA
  }

  # 9) Annotate genes (using annotateCoverageWithGenes)
  coverage_annot <- offset$coverage_data %>% mutate(end = start)
  coverage_annot <- annotateCoverageWithGenes(coverage_annot, txdb = txdb, species_annotation = org)
  if (!is.null(highlight_genes)) {
    coverage_hits <- coverage_annot %>% filter(gene_symbol %in% highlight_genes)
  } else {
    coverage_hits <- coverage_annot %>% filter(!is.na(gene_symbol)) %>% filter(abs(delta) >= gene_delta_threshold)
  }
  coverage_hits <- coverage_hits %>% group_by(gene_symbol) %>% slice_min(start, n = 1) %>% ungroup()

  # 10) Build base plot: adjust delta_weighted values to enforce min/max constraints
  if (!is.null(min_value)) {
    offset$coverage_data$delta_weighted[offset$coverage_data$delta_weighted < min_value] <- min_value
  }
  if (!is.null(max_value)) {
    offset$coverage_data$delta_weighted[offset$coverage_data$delta_weighted > max_value] <- max_value
  }
  chromosome_midpoints <- offset$coverage_data %>%
    group_by(seqnames) %>% summarize(midpoint = mean(range(cumulative_genomic_position)), .groups = 'drop')
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
      color = line_color, size = line_size, na.rm = TRUE
    ) +
    scale_color_manual(values = colors, drop = FALSE) +
    { if (nrow(data_highlight) > 0) {
      geom_rect(
        data = data_highlight,
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = VariantType),
        color = "transparent",
        alpha = variant_alpha
      )
    } else {
      NULL
    } } +
    scale_fill_manual(
      values = c("DEL" = "blue", "DUP" = "red"),
      na.value = "grey80"
    ) +
    geom_vline(
      data = offset$coverage_data %>% group_by(seqnames) %>% summarize(boundary = min(cumulative_genomic_position), .groups = "drop"),
      aes(xintercept = boundary),
      color = "grey",
      linetype = "dashed",
      alpha = 0.6
    ) +
    labs(
      title = paste0("Genome-Wide CNA ", samplename),
      x = "Chromosome",
      y = ifelse(method %in% c("loess", "fit"),
                 expression(log[2](frac(Obs, Exp)) ~ "Coverage"),
                 expression(delta ~ Coverage))
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

  # 11) Add gene labels if any gene hits are identified
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

  # 12) Optionally add G-band annotations if a band file is provided and showCNAbands is TRUE
  if (!is.null(gband_file) && showCNAbands) {
    bands_df <- fread(gband_file)
    bands_df$chr <- paste0("chr", bands_df$contig)
    band_gr <- GRanges(
      seqnames = bands_df$chr,
      ranges = IRanges(start = bands_df$start, end = bands_df$end),
      band = bands_df$name,
      giemsa = bands_df$giemsa
    )
    bands_gr <- keepStandardChromosomes(band_gr, pruning.mode = "tidy")
    var_gr <- GRanges(
      seqnames = variant_calls$seqnames,
      ranges = IRanges(start = variant_calls$start, end = variant_calls$end)
    )
    ov <- findOverlaps(var_gr, band_gr)
    if (length(ov) > 0) {
      var_band <- data.frame(
        variant_idx = queryHits(ov),
        band = mcols(band_gr)$band[subjectHits(ov)],
        giemsa = mcols(band_gr)$giemsa[subjectHits(ov)]
      )
      var_band <- var_band %>% group_by(variant_idx) %>% slice(1) %>% ungroup()
      variant_calls$band <- NA
      variant_calls$giemsa <- NA
      variant_calls$band[var_band$variant_idx] <- var_band$band
      variant_calls$giemsa[var_band$variant_idx] <- var_band$giemsa
      data_highlight$band <- variant_calls$band
      data_highlight$giemsa <- variant_calls$giemsa
      band_annotations <- data_highlight %>%
        filter(!is.na(band)) %>%
        mutate(label = band,
               x_label = (xmin + xmax) / 2)
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

#' Compute Weighted Log2 Ratios
#'
#' Computes log2 ratios of observed and reference values and calculates a weight for each observation
#' based on bin length and repeat fraction.
#'
#' @param data A \code{data.frame} containing the observed, reference, bin length, and repeat fraction columns.
#' @param observed_col Character. Name of the column with observed values.
#' @param reference_col Character. Name of the column with reference values.
#' @param bin_length_col Character. Name of the column indicating bin length.
#' @param repeat_fraction_col Character. Name of the column with repeat fraction values.
#'
#' @return The input \code{data} augmented with two new columns: \code{log2_ratio} and \code{weight}.
#'
#' @details The log2 ratio is computed as \code{log2(observed / reference)}.
#'   The weight is calculated as the product of the ratio of bin length over the median bin length and
#'   \code{(1 - repeat_fraction)}.
#' @keywords internal
#' @examples
#' df <- data.frame(observed = c(10, 20, 30),
#'                  reference = c(8, 18, 33),
#'                  bin_length = c(100, 120, 110),
#'                  repeat_fraction = c(0.1, 0.2, 0.15))
#' result <- compute_weighted_log2_ratios(df, "observed", "reference", "bin_length", "repeat_fraction")
#' print(result)
#'
#' @export
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

#' Weighted Segment Log2
#'
#' Calculates a weighted average log2 ratio for a segment given the per-observation log2 ratios and weights.
#'
#' @param segment_data A \code{data.frame} with columns \code{log2_ratio} and \code{weight}.
#'
#' @return A numeric value representing the weighted average log2 ratio. Returns \code{NA} if the total weight is zero.
#'
#' @details The function computes the sum of the weights, and if non-zero, returns the weighted average
#'   of the log2 ratios. If the total weight is zero, \code{NA} is returned.
#' @keywords internal
#' @examples
#' seg_data <- data.frame(log2_ratio = c(0.5, 1.2, -0.7),
#'                        weight = c(1, 2, 1))
#' ws <- weighted_segment_log2(seg_data)
#' print(ws)
#'
#' @export
weighted_segment_log2 <- function(segment_data) {
  total_weight <- sum(segment_data$weight, na.rm = TRUE)
  if (total_weight == 0) return(NA)
  weighted_avg <- sum(segment_data$log2_ratio * segment_data$weight, na.rm = TRUE) / total_weight
  return(weighted_avg)
}

