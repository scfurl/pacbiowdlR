
#' Plot genomic fusions using circos plots
#'
#' @description
#' Creates a circos plot of genomic fusions, with options to highlight specific fusions
#' and annotate links based on breakpoint types.  Data input should ideally be generated from the
#' 'get_fusions_from_vcf' function
#'
#' @param sv_data A data frame containing fusion data with at least the columns
#'   `chr1`, `chr2`, `start1`, and `start2`. If `annotate=TRUE`, must also include
#'   `bp1_type` and `bp2_type` columns.
#' @param highlight Integer vector of row indices to highlight with a thicker black line. Default is NULL.
#' @param annotate Logical; whether to color-code links based on fusion types (genic vs. intergenic). Default is FALSE.
#' @param title Character; title for the plot. Default is empty string.
#'
#' @return Invisibly returns NULL, but creates a circos plot as a side effect.
#'
#' @details
#' The function creates a circos plot showing genomic fusions. When `annotate=TRUE`,
#' fusions are color-coded:
#' - Red: gene-gene fusions (both breakpoints are genic)
#' - Blue: mixed fusions (one breakpoint is genic, the other intergenic)
#' - Green: intergenic fusions (both breakpoints are intergenic)
#'
#' @examples
#' \dontrun{
#' # Basic plot
#' plot_fusions(fusions)
#'
#' # Highlight the first three fusions and add annotation
#' plot_fusions(fusions, highlight = 1:3, annotate = TRUE,
#'              title = "Genomic Fusions with Highlighted Examples")
#' }
#'
#' @importFrom circlize circos.initializeWithIdeogram circos.genomicLink circos.clear
#' @importFrom graphics title legend
#'
#' @export
plot_fusions <- function(sv_data, highlight = NULL, annotate = FALSE, title = "") {
  if (!all(colnames(sv_data)[1:4] == c("chr1", "chr2", "start1", "start2"))) {
    stop("Input data must have columns 'chr1', 'chr2', 'start1', 'start2' in that order")
  }

  if (annotate && !all(c("bp1_type", "bp2_type") %in% colnames(sv_data))) {
    stop("When annotate=TRUE, input data must include 'bp1_type' and 'bp2_type' columns")
  }

  if (!is.null(highlight)) {
    highlight <- sv_data[highlight, ]
    bed3 <- highlight[, c(1, 3)]
    colnames(bed3) <- c("chr", "start")
    bed3$end <- bed3$start
    bed4 <- highlight[, c(2, 4)]
    colnames(bed4) <- c("chr", "start")
    bed4$end <- bed4$start
  }

  bed1 <- sv_data[, c(1, 3)]
  colnames(bed1) <- c("chr", "start")
  bed1$end <- bed1$start

  bed2 <- sv_data[, c(2, 4)]
  colnames(bed2) <- c("chr", "start")
  bed2$end <- bed2$start

  link_colors <- rep("grey", nrow(bed1))

  if (annotate) {
    both_genic <- sv_data$bp1_type == "genic" & sv_data$bp2_type == "genic"
    mixed_genic <- (sv_data$bp1_type == "genic" & sv_data$bp2_type != "genic") |
      (sv_data$bp1_type != "genic" & sv_data$bp2_type == "genic")
    both_intergenic <- sv_data$bp1_type != "genic" & sv_data$bp2_type != "genic"

    link_colors[both_genic] <- "#E41A1C80"
    link_colors[mixed_genic] <- "#377EB880"
    link_colors[both_intergenic] <- "#4DAF4A80"
  }

  circlize::circos.initializeWithIdeogram()
  circlize::circos.genomicLink(bed1, bed2, col = link_colors, border = NA)

  if (!is.null(highlight)) {
    circlize::circos.genomicLink(bed3, bed4, col = "black", lwd = 10)
  }

  graphics::title(title)

  if (annotate) {
    graphics::legend("bottomright",
                     legend = c("Gene-Gene Fusion", "Gene-Intergenic Fusion", "Intergenic-Intergenic"),
                     fill = c("#E41A1C80", "#377EB880", "#4DAF4A80"),
                     title = "Fusion Types", border = NA, bty = "n")
  }

  circlize::circos.clear()

  invisible(NULL)
}


#' Extract Fusions from VCF and Generate a Circos Plot
#'
#' Reads a gzipped VCF file containing structural variant data, applies several filters,
#' optionally annotates fusion breakpoints using a GTF file, and (optionally) generates a Circos plot.
#'
#' @param vcf_file Character. File path to the gzipped VCF file.
#' @param gtf_file Character or \code{NULL}. Optional file path to a GTF file for breakpoint annotation. Required when \code{annotate} is \code{TRUE}.
#' @param genome Character. Genome assembly version (default: \code{"hg38"}).
#' @param title Character. Title for the plot when \code{plot = TRUE} (default: an empty string).
#' @param thresh Numeric. Threshold for read depth filtering (default: \code{20}).
#' @param allowed Character vector. Allowed chromosomes (default: \code{paste0("chr", 1:22)}).
#' @param highlight Optional. Indices or logical vector of fusion events to highlight in the plot (default: \code{NULL}).
#' @param plot Logical. If \code{TRUE} a Circos plot will be generated (default: \code{FALSE}).
#' @param filter Character vector. Filtering criteria; supported values include \code{"pass"}, \code{"fully_spanned"}, \code{"protein_coding"}, and \code{"has_strand"}.
#' @param annotate Logical. If \code{TRUE}, fusion breakpoints will be annotated using gene coordinates (default: \code{FALSE}).
#' @param verbose Logical. if \code{TRUE}, will provide messaging.
#' @param tss_upstream Numeric. Number of base pairs upstream of the transcription start site (default: \code{2000}).
#' @param tss_downstream Numeric. Number of base pairs downstream of the transcription start site (default: \code{200}).
#'
#' @return A data frame containing fusion events with their breakpoint coordinates and (if requested) annotation details.
#'
#' @importFrom VariantAnnotation readVcf info
#' @importFrom circlize circos.initializeWithIdeogram circos.genomicLink circos.clear
#' @importFrom dplyr mutate distinct bind_rows select
#' @export

get_fusions_from_vcf <- function(vcf_file, gtf_file = NULL, genome = "hg38",
                                 title = "", thresh = 20,
                                 allowed = paste0("chr", c(1:22)),
                                 highlight = NULL, plot = FALSE,
                                 filter = c("pass", "fully_spanned", "has_strand"),
                                 annotate = FALSE, verbose = F,
                                 tss_upstream = 2000, tss_downstream = 200) {

  # If annotation is requested, preload the GTF file
  if (annotate && !is.null(gtf_file)) {
    # This will cache the GTF file if not already cached
    load_cached_gtf(gtf_file)
  }

  # Read VCF file
  if(verbose){message("Reading VCF file...")}
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

  if(verbose){message("Processing structural variants...")}
  dat <- gsub("pbsv.BND.", "", as.character(VariantAnnotation::info(vcf)$MATEID))
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
    if(verbose){message("Annotating fusion breakpoints...")}
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
    sv_data$bp1_within_exon <- bp1_annotations$within_exon
    sv_data$bp1_fiveprime_exon <- bp1_annotations$fiveprime_exon
    sv_data$bp1_threeprime_exon <- bp1_annotations$threeprime_exon
    sv_data$bp2_within_exon <- bp2_annotations$within_exon
    sv_data$bp2_fiveprime_exon <- bp2_annotations$fiveprime_exon
    sv_data$bp2_threeprime_exon <- bp2_annotations$threeprime_exon

    # For intergenic regions, add nearest genes
    sv_data$bp1_upstream_gene <- bp1_annotations$upstream_gene
    sv_data$bp1_upstream_dist <- bp1_annotations$upstream_distance
    sv_data$bp1_downstream_gene <- bp1_annotations$downstream_gene
    sv_data$bp1_downstream_dist <- bp1_annotations$downstream_distance

    sv_data$bp2_upstream_gene <- bp2_annotations$upstream_gene
    sv_data$bp2_upstream_dist <- bp2_annotations$upstream_distance
    sv_data$bp2_downstream_gene <- bp2_annotations$downstream_gene
    sv_data$bp2_downstream_dist <- bp2_annotations$downstream_distance

    if("has_strand" %in% filter){
      sv_data <- sv_data[sv_data$bp1_gene_strand %in% c("+", "-"),]
      sv_data <- sv_data[sv_data$bp2_gene_strand %in% c("+", "-"),]
    }

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
    if(verbose){message("Generating Circos plot...")}

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

#' Print Fusion List Column Names
#'
#' Diagnostic utility to display the structure of a fusion list by printing the number
#' of entries and available columns in each fusion set.
#'
#' @param fusion_list A list of data frames containing fusion data.
#' @return \code{NULL}. Prints the information directly to the console.
#' @export
#' @keywords internal
print_fusion_columns <- function(fusion_list) {
  cat("Examining fusion list structure...\n")
  for (i in 1:length(fusion_list)) {
    cat(paste0("Fusion set #", i, " has ", nrow(fusion_list[[i]]), " entries with columns:\n"))
    cat(paste(colnames(fusion_list[[i]]), collapse = ", "), "\n\n")
  }
}

#' Find Closest Related Fusions in a Reference Catalogue
#'
#' Given a set of observed (sequenced) fusions that already carry exon–level
#' breakpoint annotations, search a reference fusion database and return the
#' \code{max_results} most similar events for every observed fusion.
#'
#' Similarity is evaluated with a composite score made up of:
#' \describe{
#'   \item{Gene‐pair match}{+100 if the gene symbols match (either orientation).}
#'   \item{Strand concordance}{+50 if the strand orientations also match.}
#'   \item{Breakpoint exon proximity}{Up to +30 per partner, scaled as
#'     \eqn{30/(|d|+1)} where \eqn{d} is the absolute difference between the
#'     observed and reference exon numbers.  Both five‑prime and three‑prime
#'     nearest exons are considered.}
#'   \item{Breakpoint context}{+20 for each partner when both the observed and
#'     reference fusion breakpoints occur in exon sequence (or both in intron
#'     sequence).}
#' }
#' Only candidates with a final score \eqn{> 5} are retained.
#'
#' @param sequenced_fusions A \code{data.frame} (or \code{data.table}) produced
#'   by \code{annotate_genomic_coordinates()} containing, at minimum, the
#'   columns
#'   \itemize{
#'     \item \code{chr1}, \code{chr2}
#'     \item \code{bp1_gene}, \code{bp2_gene}
#'     \item \code{bp1_gene_strand}, \code{bp2_gene_strand}
#'     \item \code{bp1_feature},  \code{bp2_feature}
#'     \item exon‑proximity columns: \code{bp[12]\_within_exon},
#'       \code{bp[12]\_fiveprime_exon}, \code{bp[12]\_threeprime_exon}.
#'   }
#' @param db_fusions       A reference fusion catalogue as a
#'   \code{data.frame}/\code{data.table} with at least the columns
#'   \code{gene1}, \code{gene2}, \code{gene1_chro}, \code{gene2_chro},
#'   \code{gene1_strand}, \code{gene2_strand},
#'   \code{gene1_exonnumber}, \code{gene2_exonnumber}, and \code{fusion_id}.
#' @param max_results      Integer.  Maximum number of highest‑scoring database
#'   fusions returned for each sequenced fusion (default \code{5}).
#' @param verbose Logical. if \code{TRUE}, will provide messaging.
#' @return A named list.  Each element corresponds to one row in
#'   \code{sequenced_fusions} and contains up to \code{max_results} sub‑lists,
#'   each with:
#'   \itemize{
#'     \item \code{db_index} – row index in \code{db_fusions}
#'     \item \code{fusion_id} – identifier from the reference catalogue
#'     \item \code{score} – composite similarity score
#'     \item \code{gene1}, \code{gene2} – gene symbols in reference fusion
#'     \item \code{gene1_exon}, \code{gene2_exon} – reference exon numbers
#'     \item \code{match_details} – short human‑readable description
#'   }
#'
#' @details
#' Internally, two lookup hash tables are created for the reference database:
#' one keyed on gene‑pair (\code{gene1--gene2}) and one on chromosome‑pair
#' (\code{chr1_chr2}).  These drastically reduce the number of candidates that
#' must be scored for each observed fusion.
#'
#' The helper \code{calc\_score()} function (defined inside
#' \code{find_related_fusions()}) implements the scoring scheme described
#' above.  Scores \eqn{\le 5} are treated as noise and the corresponding
#' matches are discarded.
#'
#' @author Your Name
#'
#' @importFrom data.table as.data.table
#'
#' @examples
#' \donttest{
#' related <- find_related_fusions(sequenced_fusions, db_fusions, max_results = 3)
#' related[[1]]            # top matches for the first observed fusion
#' }
#' @keywords internal
#' @export

find_related_fusions <- function(sequenced_fusions, db_fusions, max_results = 5, verbose=F) {
  db_by_chr <- list()
  db_by_gene <- list()

  # Create lookup tables by chromosome pairs and gene pairs (in both orientations)
  if(verbose){message("Pre-processing database fusions...")}
  for (i in 1:nrow(db_fusions)) {
    chr1 <- db_fusions$gene1_chro[i]
    chr2 <- db_fusions$gene2_chro[i]
    gene1 <- db_fusions$gene1[i]
    gene2 <- db_fusions$gene2[i]

    # Store indices by chromosome pairs
    key_chr1 <- paste(chr1, chr2, sep = "_")
    key_chr2 <- paste(chr2, chr1, sep = "_")

    if (is.null(db_by_chr[[key_chr1]])) db_by_chr[[key_chr1]] <- integer(0)
    if (is.null(db_by_chr[[key_chr2]])) db_by_chr[[key_chr2]] <- integer(0)

    db_by_chr[[key_chr1]] <- c(db_by_chr[[key_chr1]], i)
    db_by_chr[[key_chr2]] <- c(db_by_chr[[key_chr2]], i)

    # Store indices by gene pairs
    key_gene1 <- paste(gene1, gene2, sep = "--")
    key_gene2 <- paste(gene2, gene1, sep = "--")

    if (is.null(db_by_gene[[key_gene1]])) db_by_gene[[key_gene1]] <- integer(0)
    if (is.null(db_by_gene[[key_gene2]])) db_by_gene[[key_gene2]] <- integer(0)

    db_by_gene[[key_gene1]] <- c(db_by_gene[[key_gene1]], i)
    db_by_gene[[key_gene2]] <- c(db_by_gene[[key_gene2]], i)
  }

  # Function to calculate score between a sequenced fusion and db fusion
  calc_score <- function(fusion, db_fusion) {
    score <- 0

    # Get gene names for both fusions
    bp1_gene <- if (!is.na(fusion$bp1_gene)) fusion$bp1_gene else NULL
    bp2_gene <- if (!is.na(fusion$bp2_gene)) fusion$bp2_gene else NULL

    # Score based on gene name matches (highest priority)
    if (!is.null(bp1_gene) && !is.null(bp2_gene)) {
      if ((bp1_gene == db_fusion$gene1 && bp2_gene == db_fusion$gene2) ||
          (bp1_gene == db_fusion$gene2 && bp2_gene == db_fusion$gene1)) {
        score <- score + 100
      } else {
        # If genes don't match at all, give a very low base score
        # This should prevent high scores for non-matching genes
        score <- score + 1
      }
    }

    # Score based on strand matches (only if gene names match)
    if (score > 50) {  # Only count strand if genes match
      bp1_strand <- if (!is.na(fusion$bp1_gene_strand)) fusion$bp1_gene_strand else NULL
      bp2_strand <- if (!is.na(fusion$bp2_gene_strand)) fusion$bp2_gene_strand else NULL

      if (!is.null(bp1_strand) && !is.null(bp2_strand)) {
        if ((bp1_strand == db_fusion$gene1_strand && bp2_strand == db_fusion$gene2_strand) ||
            (bp1_strand == db_fusion$gene2_strand && bp2_strand == db_fusion$gene1_strand)) {
          score <- score + 50
        }
      }
    }

    # Only calculate exon scores if gene names match
    if (score > 50) {
      # Score based on exon proximity
      bp1_exon <- NULL
      if (!is.na(fusion$bp1_fiveprime_exon) && !is.na(fusion$bp1_threeprime_exon)) {
        bp1_exon <- c(as.numeric(fusion$bp1_fiveprime_exon), as.numeric(fusion$bp1_threeprime_exon))
      }

      bp2_exon <- NULL
      if (!is.na(fusion$bp2_fiveprime_exon) && !is.na(fusion$bp2_threeprime_exon)) {
        bp2_exon <- c(as.numeric(fusion$bp2_fiveprime_exon), as.numeric(fusion$bp2_threeprime_exon))
      }

      if (!is.null(bp1_exon) && !is.na(db_fusion$gene1_exonnumber)) {
        exon_dist <- min(abs(bp1_exon[1] - db_fusion$gene1_exonnumber),
                         abs(bp1_exon[2] - db_fusion$gene1_exonnumber))
        score <- score + 30 / (exon_dist + 1)
      }

      if (!is.null(bp2_exon) && !is.na(db_fusion$gene2_exonnumber)) {
        exon_dist <- min(abs(bp2_exon[1] - db_fusion$gene2_exonnumber),
                         abs(bp2_exon[2] - db_fusion$gene2_exonnumber))
        score <- score + 30 / (exon_dist + 1)
      }

      # Consider breakpoint types (exon vs. intron)
      bp1_feature <- fusion$bp1_feature
      bp2_feature <- fusion$bp2_feature

      if (!is.na(bp1_feature) && !is.na(bp2_feature)) {
        bp1_is_exon <- grepl("exon", bp1_feature)
        bp2_is_exon <- grepl("exon", bp2_feature)

        # Check if db fusion has breakpoints in exons
        db_bp1_is_exon <- grepl("exon", db_fusion$gene1_exon)
        db_bp2_is_exon <- grepl("exon", db_fusion$gene2_exon)

        # Add score for matching breakpoint types
        if (bp1_is_exon == db_bp1_is_exon) score <- score + 20
        if (bp2_is_exon == db_bp2_is_exon) score <- score + 20
      }
    }

    return(score)
  }

  # Initialize results list
  results <- vector("list", nrow(sequenced_fusions))
  names(results) <- as.character(1:nrow(sequenced_fusions))

  # Process sequenced fusions
  if(verbose){message("Processing sequenced fusions...")}
  for (idx in 1:nrow(sequenced_fusions)) {
    fusion <- sequenced_fusions[idx, ]

    # Extract information
    chr1 <- fusion$chr1
    chr2 <- fusion$chr2
    bp1_gene <- if (!is.na(fusion$bp1_gene)) fusion$bp1_gene else NULL
    bp2_gene <- if (!is.na(fusion$bp2_gene)) fusion$bp2_gene else NULL

    # Match by gene names first (most specific)
    match_indices <- integer(0)
    if (!is.null(bp1_gene) && !is.null(bp2_gene)) {
      gene_key <- paste(bp1_gene, bp2_gene, sep = "--")
      gene_matches <- db_by_gene[[gene_key]]
      if (!is.null(gene_matches) && length(gene_matches) > 0) {
        match_indices <- gene_matches
      }
    }

    # If no gene matches, try chromosome matches
    if (length(match_indices) == 0) {
      chr_key <- paste(chr1, chr2, sep = "_")
      chr_matches <- db_by_chr[[chr_key]]
      if (!is.null(chr_matches) && length(chr_matches) > 0) {
        match_indices <- chr_matches
      }
    }

    # If still no matches, skip this fusion
    if (length(match_indices) == 0) {
      results[[as.character(idx)]] <- list()
      next
    }

    # Calculate scores for matching fusions
    scores <- numeric(length(match_indices))
    for (j in seq_along(match_indices)) {
      db_idx <- match_indices[j]
      scores[j] <- calc_score(fusion, db_fusions[db_idx, ])
    }

    # Find top matches
    if (length(scores) == 0) {
      results[[as.character(idx)]] <- list()
      next
    }

    # Sort matches by score and get top results
    ranked_order <- order(scores, decreasing = TRUE)
    top_indices <- ranked_order[1:min(max_results, length(ranked_order))]

    # Store results for this fusion
    fusion_results <- list()
    for (i in seq_along(top_indices)) {
      match_idx <- match_indices[top_indices[i]]
      match_score <- scores[top_indices[i]]

      if (match_score > 5) {  # Only include matches with meaningful scores
        match_fusion <- db_fusions[match_idx, ]
        fusion_results[[i]] <- list(
          db_index = match_idx,
          fusion_id = match_fusion$fusion_id,
          score = match_score,
          gene1 = match_fusion$gene1,
          gene2 = match_fusion$gene2,
          gene1_exon = match_fusion$gene1_exonnumber,
          gene2_exon = match_fusion$gene2_exonnumber,
          match_details = sprintf(
            "Match for fusion %d: %s (score: %.2f) - Gene1: %s (exon %s), Gene2: %s (exon %s)",
            idx, match_fusion$fusion_id, match_score,
            match_fusion$gene1, match_fusion$gene1_exonnumber,
            match_fusion$gene2, match_fusion$gene2_exonnumber
          )
        )
      }
    }

    results[[as.character(idx)]] <- fusion_results
  }

  return(results)
}


#' Print fusion matches in a readable format
#'
#' @description
#' Formats and displays the results from the fusion matching function in a human-readable
#' summary format, showing all matches found for each fusion.
#'
#' @param results A list of fusion match results, as returned by \code{find_related_fusions()}.
#'   Each element corresponds to a sequenced fusion and contains nested lists of matched
#'   database fusions with scores and details.
#'
#' @return Invisibly returns NULL. As a side effect, prints formatted match results
#'   to the console.
#'
#' @details
#' For each sequenced fusion, this function prints:
#' \itemize{
#'   \item A header showing the fusion index and match count
#'   \item No matches found message (if applicable)
#'   \item For each match: index, fusion ID, score, and gene/exon details
#' }
#'
#' The function handles empty results and properly formats each match's details
#' based on the structure provided by \code{find_related_fusions()}.
#'
#' @examples
#' \dontrun{
#' # Find related fusions
#' results <- find_related_fusions(sequenced_fusions, db_fusions)
#'
#' # Print matches in readable format
#' print_fusion_matches(results)
#'
#' # Print matches for just one fusion
#' print_fusion_matches(results[1])
#' }
#'
#' @seealso \code{\link{find_related_fusions}} for generating the match results
#'
#' @export
print_fusion_matches <- function(results) {
  cat("Fusion Matches Summary:\n")
  cat("=======================\n\n")

  for (fusion_idx in names(results)) {
    matches <- results[[fusion_idx]]

    if (length(matches) == 0) {
      cat(sprintf("Fusion %s: No matches found\n\n", fusion_idx))
      next
    }

    cat(sprintf("Fusion %s: Found %d matches\n", fusion_idx, length(matches)))

    for (i in 1:length(matches)) {
      match <- matches[[i]]
      cat(sprintf("  %d. %s\n", i, match$match_details))
    }
    cat("\n")
  }
}

#' Get the NeoSplice Fusion Database
#'
#' Retrieves (and caches) the NeoSplice fusion database. If a file path is provided, it is loaded and parsed;
#' otherwise, a default database (bundled with the package) is used.
#'
#' @param db_file Character or \code{NULL}. Optional file path to the NeoSplice database file.
#' @param force_reload Logical. If \code{TRUE}, forces reloading of the database even if cached (default: \code{FALSE}).
#' @param verbose Logical. if \code{TRUE}, will provide messaging.
#' @return A data frame containing the NeoSplice fusion database entries.
#'
#' @importFrom digest digest
#' @keywords internal
get_neosplice_db <- function(db_file = system.file("extdata", "fusions/NeoSplice_hg38_inframe_fusion.txt.gz",
                                                   package = "pacbiowdlR"), force_reload = FALSE, verbose=F) {
  db_key <- "neosplice_db"

  # Use custom database if provided
  if (!is.null(db_file)) {
    if (!file.exists(db_file)) {
      stop("Database file not found: ", db_file)
    }
    db_key <- paste0("neosplice_db_", digest::digest(db_file))
  } else {
    # Use default database from package resources
    db_file <- system.file("extdata", "fusions", "NeoSplice_hg38_inframe_fusion.txt.gz",
                           package = "yourpackage")  # Replace "yourpackage" with your package name
  }

  # Check if already cached
  if (!force_reload && exists(db_key, envir = .fusion_pkg_env)) {
    return(get(db_key, envir = .fusion_pkg_env))
  }

  if(verbose){message("Loading NeoSplice database...")}
  db <- parse_neosplice_database(db_file)

  # Cache the database
  assign(db_key, db, envir = .fusion_pkg_env)

  return(db)
}


#' Parse NeoSplice fusion database
#'
#' Specifically designed parser for the NeoSplice fusion database format.
#' This function correctly handles the column naming and format specific to NeoSplice.
#'
#' @param file_path Path to the fusion database file (can be gzipped, with .gz extension)
#' @param verbose Logical. if \code{TRUE}, will provide messaging.
#' @return A data frame with standardized fusion database information
#' @export
#' @examples
#' # Parse a NeoSplice database
#' db_fusions <- parse_neosplice_database("path/to/NeoSplice_hg38_inframe_fusion.txt.gz")
#' @keywords internal
#' # Examine the first few rows
#' head(db_fusions)
parse_neosplice_database <- function(file_path, verbose=F) {
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
    if(verbose){message("Error parsing NeoSplice database: ", e$message)}

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

    if(verbose){message("Header: ", header)}
    if(verbose){message("Sample lines: ")}
    for (line in sample_lines) {
      if(verbose){message(line)}
    }

    stop("Failed to parse NeoSplice database", call. = FALSE)
  })
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
#' @param verbose Logical. if \code{TRUE}, will provide messaging.
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
#' @keywords internal
#' @export
get_fusion_sequences <- function(fusion_df, n_bp = 100, verbose=F) {
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
      if(verbose){warning(paste("Row", i, "has NA in bp1_gene_strand. Defaulting to '*' strand."))}
      strand1 <- "*"
    }

    if (strand1 == "-") {
      # For + strand, take downstream of breakpoint
      range1 <- GRanges(seqnames = chr1,
                        ranges = IRanges(start = start1, end = start1 + n_bp - 1),
                        strand = "+")
      seq1 <- as.character(getSeq(genome, range1))
    } else {
      # For - strand, take upstream of breakpoint, then reverse complement
      range1 <- GRanges(seqnames = chr1,
                        ranges = IRanges(start = start1 - n_bp + 1, end = start1),
                        strand = "+")
      seq1 <- as.character(getSeq(genome, range1))
    }

    # Handle second breakpoint strand (default to "+" if NA)
    if (is.na(strand2)) {
      if(verbose){warning(paste("Row", i, "has NA in bp2_gene_strand. Defaulting to '*' strand."))}
      strand2 <- "*"
    }

    if (strand2 == "-") {
      # For + strand, take downstream of breakpoint
      range2 <- GRanges(seqnames = chr2,
                        ranges = IRanges(start = start2, end = start2 + n_bp - 1),
                        strand = "-")
      seq2 <- as.character(getSeq(genome, range2))
    } else {
      # For - strand, take upstream of breakpoint, then reverse complement
      range2 <- GRanges(seqnames = chr2,
                        ranges = IRanges(start = start2 - n_bp + 1, end = start2),
                        strand = "+")
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
#' @param verbose Logical. if \code{TRUE}, will provide messaging.
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
#' @keywords internal
#' @export
create_fusion_fasta <- function(fusion_df, file_path = "fusion_sequences.fasta", include_genes = TRUE, verbose=F) {
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
  if(verbose){message(paste("Created FASTA file:", file_path))}
}


