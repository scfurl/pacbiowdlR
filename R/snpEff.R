library(tidyverse)
library(vroom)
library(rsnps)
library(gt)
# remotes::install_github("ropensci/rsnps")
library(future.apply)  # For parallelization
library(tidyr)

# Plan for parallel execution (adjust workers as needed)


plot_snpEff <- function(file, cores = parallel::detectCores() - 1){
  options(future.globals.maxSize = 8 * 1024^3)  # 2 GiB (adjust as needed)

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

  plan(multisession, workers = cores)
  ann <- vroom(file, guess_max = 1000000)

  ann_annotations <- ann$ANN

  # table based on SnpEFF documentation
  ann_types <- data.frame(impact = c("HIGH","HIGH","HIGH","HIGH","HIGH","HIGH","HIGH","HIGH","HIGH","HIGH","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE","MODERATE","LOW","LOW","LOW","LOW","LOW","LOW","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER","MODIFIER"),
                          ontology_term = c("chromosome_number_variation","exon_loss_variant","frameshift_variant","rare_amino_acid_variant","splice_acceptor_variant","splice_donor_variant","start_lost","stop_gained","stop_lost","transcript_ablation","3_prime_UTR_truncation&exon_loss","5_prime_UTR_truncation&exon_loss_variant","coding_sequence_variant-moderate","conservative_inframe_deletion","conservative_inframe_insertion","disruptive_inframe_deletion","disruptive_inframe_insertion","missense_variant","regulatory_region_ablation","splice_region_variant-moderate","TFBS_ablation","5_prime_UTR_premature_start_codon_gain_variant","initiator_codon_variant","splice_region_variant-low","start_retained","stop_retained_variant","synonymous_variant","3_prime_UTR_variant","5_prime_UTR_variant","coding_sequence_variant-modifier","conserved_intergenic_variant","conserved_intron_variant","downstream_gene_variant","exon_variant","feature_elongation","feature_truncation","gene_variant","intergenic_region","intragenic_variant","intron_variant","mature_miRNA_variant","miRNA","NMD_transcript_variant","non_coding_transcript_exon_variant","non_coding_transcript_variant","regulatory_region_amplification","regulatory_region_variant","TF_binding_site_variant","TFBS_amplification","transcript_amplification","transcript_variant","upstream_gene_variant"))


  # **Parallelized variant counting using `future_lapply()`**
  ontology_terms <- ann_types$ontology_term  # Avoid passing full tibble
  n_variants <- future_lapply(ontology_terms, function(term) {
    sum(str_detect(ann_annotations, term), na.rm = TRUE)
  })

  # Convert result back to tibble
  ann_types$n_variants <- unlist(n_variants)

  # Sort and transform data
  ann_types <- ann_types %>%
    arrange(n_variants) %>%
    mutate(
      ontology_term = factor(ontology_term, levels = ontology_term),
      impact = factor(impact, levels = c("HIGH","MODERATE","LOW","MODIFIER")),
      log_n_variants = log10(n_variants + 1)  # Avoid log(0)
    ) %>%
    filter(n_variants > 0)



  p1 <- ann_types|>
    ggplot(aes(ontology_term, log_n_variants, fill = impact)) +
    geom_col() +
    coord_flip() +
    theme_minimal() +
    theme(axis.title.x=element_text(size=14),
          axis.title.y=element_text(size=12),
          axis.title=element_text(size=14),
          axis.text=element_text(size=12),
          plot.title=element_text(size=18))+
    ylab("\nlog(number of annotated variants)")+ xlab("")+ ylim(0,8)+
    ggtitle("Number of variants showing each snpEFF annotation (log scale)")+
    geom_text(aes(label=n_variants), position=position_stack(), hjust=-0.5)+
    scale_fill_manual(values = c(red,yel,blu,gre))
  p1
}
