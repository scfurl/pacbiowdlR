# Load required libraries
library(circlize)
library(VariantAnnotation)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(IRanges)
library(GenomeInfoDb)

# Function to read gzipped VCF and generate a Circos plot with annotations
plot_circos_from_vcf <- function(vcf_file, gtf_file = NULL, genome = "hg38",
                                 title = "", thresh = 20,
                                 allowed = paste0("chr", c(1:22)),
                                 highlight = NULL, return_data = TRUE,
                                 annotate = FALSE,
                                 tss_upstream = 2000, tss_downstream = 200) {
  # Check if file exists
  if (!file.exists(vcf_file)) {
    stop(paste("Error: File not found:", vcf_file))
  }

  # Check if GTF file exists if annotation is requested
  if (annotate && is.null(gtf_file)) {
    stop("Error: GTF file required for annotation. Please provide gtf_file parameter.")
  }

  if (annotate && !file.exists(gtf_file)) {
    stop(paste("Error: GTF file not found:", gtf_file))
  }

  # Read VCF file
  message("Reading VCF file...")
  vcf <- readVcf(vcf_file, genome = genome)
  vcf <- vcf[vcf@info$SVTYPE=="BND",]
  vcf <- vcf[vcf@fixed$FILTER=="PASS"]
  vcf <- vcf[vcf@info$NotFullySpanned=="FALSE"]
  vcf <- vcf[grepl("protein_coding", sapply(vcf@info$BCSQ, function(bs){ifelse(length(bs)==0, NA, bs[1])})),]
  vcf <- vcf[vcf@assays@data$DP[,1]>thresh,]

  message("Processing structural variants...")
  dat <- gsub("pbsv.BND.", "", as.character(info(vcf)$MATEID))
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
    message("Annotating fusion breakpoints...")
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

    # For intergenic regions, add nearest genes
    sv_data$bp1_upstream_gene <- bp1_annotations$upstream_gene
    sv_data$bp1_upstream_dist <- bp1_annotations$upstream_distance
    sv_data$bp1_downstream_gene <- bp1_annotations$downstream_gene
    sv_data$bp1_downstream_dist <- bp1_annotations$downstream_distance

    sv_data$bp2_upstream_gene <- bp2_annotations$upstream_gene
    sv_data$bp2_upstream_dist <- bp2_annotations$upstream_distance
    sv_data$bp2_downstream_gene <- bp2_annotations$downstream_gene
    sv_data$bp2_downstream_dist <- bp2_annotations$downstream_distance

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
  }

  if (return_data) {
    return(sv_data)
  } else {
    message("Generating Circos plot...")

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
}
#
# # This is the annotate_genomic_coordinates function required for annotation
# # Copy from previous code - this is a placeholder for the full function
# annotate_genomic_coordinates <- function(coordinates, genome, gtffile,
#                                          tss_upstream = 2000, tss_downstream = 200,
#                                          gene_types = NULL, transcript_types = NULL,
#                                          return_all_as_df = TRUE) {
#   # Check if required packages are installed
#   required_packages <- c("rtracklayer", "GenomicRanges", "IRanges", "GenomeInfoDb")
#   missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
#
#   if (length(missing_packages) > 0) {
#     stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
#          ". Please install with BiocManager::install()")
#   }
#
#   # Load necessary libraries
#   suppressPackageStartupMessages({
#     library(rtracklayer)
#     library(GenomicRanges)
#     library(IRanges)
#     library(GenomeInfoDb)
#   })
#
#   # Check if gtffile exists
#   if (!file.exists(gtffile)) {
#     stop("GTF file not found: ", gtffile)
#   }
#
#   # Process input coordinates into a standard data frame
#   if (is.data.frame(coordinates)) {
#     # Check if required columns exist
#     if (!all(c("chr", "pos") %in% colnames(coordinates))) {
#       stop("Input data frame must contain 'chr' and 'pos' columns")
#     }
#     coord_df <- coordinates[, c("chr", "pos")]
#   } else if (is.list(coordinates)) {
#     if (all(c("chr", "pos") %in% names(coordinates))) {
#       # List with chr and pos vectors
#       if (length(coordinates$chr) != length(coordinates$pos)) {
#         stop("'chr' and 'pos' vectors must have the same length")
#       }
#       coord_df <- data.frame(chr = coordinates$chr, pos = coordinates$pos)
#     } else {
#       # List of individual coordinates
#       coord_df <- do.call(rbind, lapply(coordinates, function(coord) {
#         if (is.list(coord) && all(c("chr", "pos") %in% names(coord))) {
#           data.frame(chr = coord$chr, pos = coord$pos)
#         } else {
#           stop("Each list item must contain 'chr' and 'pos' elements")
#         }
#       }))
#     }
#   } else {
#     # Single coordinate pair
#     if (is.character(coordinates) && is.numeric(pos)) {
#       coord_df <- data.frame(chr = coordinates, pos = pos)
#     } else {
#       stop("Coordinates must be provided as a data frame, list, or chr/pos pair")
#     }
#   }
#
#   # Add row identifier if there isn't one
#   if (!"id" %in% colnames(coord_df)) {
#     coord_df$id <- seq_len(nrow(coord_df))
#   }
#
#   # Import GTF file
#   message("Importing GTF file...")
#   gtf_data <- import(gtffile)
#
#   # Extract genes from GTF
#   genes <- gtf_data[gtf_data$type == "gene"]
#
#   # Ensure genes have gene_name or gene_id
#   if ("gene_name" %in% names(mcols(genes))) {
#     genes$symbol <- genes$gene_name
#   } else if ("gene_id" %in% names(mcols(genes))) {
#     genes$symbol <- genes$gene_id
#   } else {
#     genes$symbol <- rep(NA, length(genes))
#   }
#
#   # Create dummy annotation dataframe for demonstration
#   # In a real implementation, this would use the full annotation logic
#   result_df <- data.frame(
#     id = coord_df$id,
#     chr = coord_df$chr,
#     position = coord_df$pos,
#     genome = rep(genome, nrow(coord_df)),
#     location_type = sample(c("genic", "intergenic"), nrow(coord_df), replace = TRUE),
#     feature = sample(c("exon", "intron", "5' UTR", "3' UTR"), nrow(coord_df), replace = TRUE),
#     gene_symbol = sample(genes$symbol[1:100], nrow(coord_df), replace = TRUE),
#     gene_strand = sample(c("+", "-"), nrow(coord_df), replace = TRUE),
#     gene_type = sample(c("protein_coding", "lincRNA"), nrow(coord_df), replace = TRUE),
#     upstream_gene = sample(genes$symbol[1:100], nrow(coord_df), replace = TRUE),
#     upstream_distance = sample(1000:50000, nrow(coord_df), replace = TRUE),
#     downstream_gene = sample(genes$symbol[1:100], nrow(coord_df), replace = TRUE),
#     downstream_distance = sample(1000:50000, nrow(coord_df), replace = TRUE)
#   )
#
#   message("Annotation complete")
#   return(result_df)
# }
#
# # Example usage:
# # 1. To get annotated fusion data:
# # fusions <- plot_circos_from_vcf("fusions.vcf.gz", "annotations.gtf", annotate = TRUE, return_data = TRUE)
# #
# # 2. To plot with annotations:
# # plot_circos_from_vcf("fusions.vcf.gz", "annotations.gtf", annotate = TRUE,
# #                    highlight = c(1, 3), return_data = FALSE,
# #                    title = "Fusion Breakpoints with Annotation")
# # # Load required libraries
# # library(circlize)
# # library(VariantAnnotation)
# # library(dplyr)
# #
# # # Function to read gzipped VCF and generate a Circos plot
# # plot_circos_from_vcf <- function(vcf_file,title = "", thresh = 20, allowed = paste0("chr", c(1:22)), highlight = NULL, return_data = T, annotate = F){
# #
# #   # Check if file exists
# #   if (!file.exists(vcf_file)) {
# #     stop(paste("Error: File not found:", vcf_file))
# #   }
# #
# #   # Read VCF file
# #   vcf <- readVcf(vcf_file, genome = "GRCh38")
# #   vcf <- vcf[vcf@info$SVTYPE=="BND",]
# #   vcf <- vcf[vcf@fixed$FILTER=="PASS"]
# #   vcf <- vcf[vcf@info$NotFullySpanned=="FALSE"]
# #   vcf <- vcf[grepl("protein_coding", sapply(vcf@info$BCSQ, function(bs){ifelse(length(bs)==0, NA, bs[1])})),]
# #   vcf <- vcf[vcf@assays@data$DP[,1]>thresh,]
# #
# #   dat <- gsub("pbsv.BND.", "", as.character(info(vcf)$MATEID))
# #   parts <- strsplit(dat, "-")
# #   sv_data <- do.call(rbind, parts)
# #   sv_data <-data.frame(chr1=strsplit(sv_data[,1], ":") %>% sapply("[[", 1),
# #              chr2=strsplit(sv_data[,2], ":") %>% sapply("[[", 1),
# #              start1=strsplit(sv_data[,1], ":") %>% sapply("[[", 2),
# #              start2=strsplit(sv_data[,2], ":") %>% sapply("[[", 2))
# #   sv_data <- sv_data[sv_data$chr1 %in% allowed & sv_data$chr2 %in% allowed, ]
# #
# #   sv_data <- sv_data %>%
# #     mutate(pair_key = pmin(chr1, chr2), pair_val = pmax(chr1, chr2),
# #            start_key = pmin(start1, start2), start_val = pmax(start1, start2)) %>%
# #     distinct(pair_key, pair_val, start_key, start_val, .keep_all = TRUE) %>%
# #     dplyr::select(-pair_key, -pair_val, -start_key, -start_val)
# #
# #
# #   # Remove NA values (only keep inter-chromosomal events)
# #   sv_data <- sv_data[!is.na(sv_data$chr2) & !is.na(sv_data$start2), ]
# #
# #
# #   if(return_data){
# #     return(sv_data)
# #   } else {
# #     highlight <- sv_data[highlight, ]
# #     bed3<-highlight[,c(1,3)]
# #     colnames(bed3) <- c("chr", "start")
# #     bed3$start = as.numeric(bed3$start )
# #     bed3$end <- bed3$start
# #     bed4 = as.data.frame(highlight[,c(2,4)])
# #     colnames(bed4) <- c("chr", "start")
# #     bed4$start = as.numeric(bed4$start )
# #     bed4$end <- bed4$start
# #
# #     bed1 = sv_data[,c(1,3)]
# #     colnames(bed1) <- c("chr", "start")
# #     bed1$start = as.numeric(bed1$start )
# #     bed1$end <- bed1$start
# #     bed2 = as.data.frame(sv_data[,c(2,4)])
# #     colnames(bed2) <- c("chr", "start")
# #     bed2$start = as.numeric(bed2$start )
# #     bed2$end <- bed2$start
# #
# #
# #     circos.initializeWithIdeogram()
# #     circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5),
# #                        border = NA)
# #     circos.genomicLink(bed3, bed4, col = "black",
# #                        lwd = 10)
# #     circos.clear()
# #     title(title)
# #   }
# #
# #
# # }
# #
