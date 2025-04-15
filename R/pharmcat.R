#' @export
#' @keywords internal
pharmcat_json <- function(json_file, genes = c("CYP2B6", "CYP2C9",   "CYP3A5",    "G6PD",   "TPMT", "NUDT15")){

  data <- fromJSON(json_file, flatten = TRUE)

  # Extract relevant gene information
  genes_data <- data$genes$CPIC  # Adjust the path if needed
  to_eval <- names(genes_data)[names(genes_data) %in% genes]
  # Convert to a dataframe
  gene_list <- lapply(to_eval, function(gene) {
    gene_info <- genes_data[[gene]]
    message(paste0("evaluating gene: ", gene))
    # Extract relevant fields
    data.frame(
      Gene = gene_info$geneSymbol,
      Phenotype = paste(gene_info$sourceDiplotypes$phenotypes, collapse = "; "),
      Allele1 = gene_info$sourceDiplotypes$allele1.name,
      Function1 = gene_info$sourceDiplotypes$allele1.function,
      Allele2 = gene_info$sourceDiplotypes$allele2.name,
      Function2 = gene_info$sourceDiplotypes$allele2.function,
      RelatedDrugs = paste(gene_info$relatedDrugs$name, collapse = ", "),
      stringsAsFactors = FALSE
    )
  })

  # Combine into a single dataframe
  df_genes <- bind_rows(gene_list)

  df_genes
}
