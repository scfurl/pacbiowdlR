#' Plot SnpEff Annotation Impact Distribution
#'
#' @description
#' Creates a horizontal bar chart visualizing the distribution of SnpEff variant annotations,
#' categorized by impact level (HIGH, MODERATE, LOW, MODIFIER). The chart displays the number
#' of variants for each annotation type on a logarithmic scale.
#'
#' @param file Character; path to a tab-delimited SnpEff annotation file.
#' @param cores Integer; number of CPU cores to use for parallel processing (default: one less than total cores).
#'
#' @return A ggplot2 object with a horizontal bar chart of SnpEff annotations.
#'
#' @details
#' This function parses a SnpEff annotation file and counts occurrences of each annotation type.
#' Annotations are categorized by impact level and displayed on a logarithmic scale. The function
#' uses parallel processing to improve performance when handling large annotation files.
#'
#' The color scheme is based on Wong's colorblind-friendly palette (Wong, B. Points of view:
#' Color blindness. Nat Methods, 2011).
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' plot_snpEff("path/to/snpeff_annotations.txt.gz")
#'
#' # Specify number of cores for parallel processing
#' plot_snpEff("path/to/snpeff_annotations.txt.gz", cores = 4)
#' }
#'
#' @importFrom vroom vroom
#' @importFrom future.apply future_lapply
#' @importFrom future plan multisession
#' @importFrom parallel detectCores
#' @importFrom stringr str_detect
#' @importFrom dplyr arrange mutate filter
#' @importFrom ggplot2 ggplot aes geom_col coord_flip theme_minimal theme
#' @importFrom ggplot2 element_text ylab xlab ylim ggtitle geom_text
#' @importFrom ggplot2 scale_fill_manual
#' @importFrom tidyr %>%
#'
#' @export
plot_snpEff <- function(file, cores = parallel::detectCores() - 1) {
  options(future.globals.maxSize = 32 * 1024^3)  # 8 GiB (adjust as needed)

  # Wong, B. Points of view: Color blindness. Nat Methods (2011).
  bla <- '#000000'
  blu <- '#0072b2'
  grb <- '#56b4e9'
  lir <- '#cc79a7'
  gre <- '#009e73'
  red <- '#d55e00'
  org <- '#e69f00'
  yel <- '#f0e442'
  gry <- '#BBBBBB'

  future::plan(future::multisession, workers = cores)

  ann <- vroom::vroom(file, guess_max = 2000, comment = "##", show_col_types = FALSE)
  if("ANN" %in% colnames(ann)){
    ann_annotations <- ann$ANN
  } else if("INFO" %in% colnames(ann)) {
    ann_annotations <- strsplit(ann$INFO, ";") %>% sapply("[[", 3) %>% strsplit(., "\\|") %>% sapply("[[", 2)
  }

  # Table based on SnpEFF documentation
  ann_types <- data.frame(
    impact = c("HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH", "HIGH",
               "MODERATE", "MODERATE", "MODERATE", "MODERATE", "MODERATE", "MODERATE", "MODERATE",
               "MODERATE", "MODERATE", "MODERATE", "MODERATE", "LOW", "LOW", "LOW", "LOW", "LOW",
               "LOW", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER",
               "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER",
               "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER",
               "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER", "MODIFIER"),
    ontology_term = c("chromosome_number_variation", "exon_loss_variant", "frameshift_variant",
                      "rare_amino_acid_variant", "splice_acceptor_variant", "splice_donor_variant",
                      "start_lost", "stop_gained", "stop_lost", "transcript_ablation",
                      "3_prime_UTR_truncation&exon_loss", "5_prime_UTR_truncation&exon_loss_variant",
                      "coding_sequence_variant-moderate", "conservative_inframe_deletion",
                      "conservative_inframe_insertion", "disruptive_inframe_deletion",
                      "disruptive_inframe_insertion", "missense_variant", "regulatory_region_ablation",
                      "splice_region_variant-moderate", "TFBS_ablation",
                      "5_prime_UTR_premature_start_codon_gain_variant", "initiator_codon_variant",
                      "splice_region_variant-low", "start_retained", "stop_retained_variant",
                      "synonymous_variant", "3_prime_UTR_variant", "5_prime_UTR_variant",
                      "coding_sequence_variant-modifier", "conserved_intergenic_variant",
                      "conserved_intron_variant", "downstream_gene_variant", "exon_variant",
                      "feature_elongation", "feature_truncation", "gene_variant", "intergenic_region",
                      "intragenic_variant", "intron_variant", "mature_miRNA_variant", "miRNA",
                      "NMD_transcript_variant", "non_coding_transcript_exon_variant",
                      "non_coding_transcript_variant", "regulatory_region_amplification",
                      "regulatory_region_variant", "TF_binding_site_variant", "TFBS_amplification",
                      "transcript_amplification", "transcript_variant", "upstream_gene_variant")
  )

  # Parallelized variant counting using future_lapply()
  ontology_terms <- ann_types$ontology_term  # Avoid passing full tibble
  n_variants <- future.apply::future_lapply(ontology_terms, function(term) {
    sum(stringr::str_detect(ann_annotations, term), na.rm = TRUE)
  })

  # Convert result back to data frame
  ann_types$n_variants <- unlist(n_variants)

  # Sort and transform data
  ann_types <- ann_types %>%
    dplyr::arrange(n_variants) %>%
    dplyr::mutate(
      ontology_term = factor(ontology_term, levels = ontology_term),
      impact = factor(impact, levels = c("HIGH", "MODERATE", "LOW", "MODIFIER")),
      log_n_variants = log10(n_variants + 1)  # Avoid log(0)
    ) %>%
    dplyr::filter(n_variants > 0)

  # Create the plot
  p1 <- ann_types %>%
    ggplot2::ggplot(ggplot2::aes(ontology_term, log_n_variants, fill = impact)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_text(size = 10),
      axis.title.y = ggplot2::element_text(size = 9),
      axis.title = ggplot2::element_text(size = 10),
      axis.text = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(size = 12)
    ) +
    ggplot2::ylab("\nlog(number of annotated variants)") +
    ggplot2::xlab("") +
    ggplot2::ylim(0, 8) +
    ggplot2::ggtitle("snpEFF annotation of variants") +
    ggplot2::geom_text(ggplot2::aes(label = n_variants),
                       position = ggplot2::position_stack(), hjust = -0.5) +
    ggplot2::scale_fill_manual(values = c(red, yel, blu, gre))

  return(p1)
}
