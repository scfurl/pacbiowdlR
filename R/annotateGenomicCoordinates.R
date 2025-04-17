#' Annotate Multiple Genomic Coordinates Using GTF
#'
#' @description
#' Improved function to annotate genomic coordinates based on a GTF file,
#' with optional canonical‑transcript exon numbering information.
#'
#' **Key speed‑ups**
#' * Canonical transcript and exon information is now pre‑computed once in `preload_gtf()` and cached.
#' * Robust checks ensure a non‑NULL `canonical_exons_by_gene` object; if unavailable, the cache is rebuilt.
#'
#' @param coordinates Data frame with columns `chr` and `pos`.
#' @param genome Character; genome build name (e.g., `"hg38"`, `"hg19"`, `"mm10"`).
#' @param gtffile Character; path to the GTF annotation file.
#' @param tss_upstream Numeric; bases upstream of TSS to define as promoter (default: 2000).
#' @param tss_downstream Numeric; bases downstream of TSS to include in promoter (default: 200).
#' @param cache_gtf Logical; whether to cache the GTF data between calls (default: TRUE).
#' @param include_exon_info Logical; whether to compute canonical‑transcript exon numbers (default: FALSE).
#'
#' @return A data frame with annotation results for each coordinate. If `include_exon_info = TRUE`,
#' extra columns `within_exon`, `fiveprime_exon`, and `threeprime_exon` are included.
#' @export
#' @keywords internal
#'
#' @importFrom GenomicRanges GRanges seqnames findOverlaps strand width
#' @importFrom IRanges IRanges start end
#' @importFrom rtracklayer import
#' @importFrom GenomeInfoDb seqlevels
#' @importFrom S4Vectors mcols
#' @importFrom digest digest
annotate_genomic_coordinates <- function(
    coordinates, genome, gtffile,
    tss_upstream = 2000, tss_downstream = 200,
    cache_gtf = TRUE,
    include_exon_info = TRUE
) {
  ## ---------------- Validations ---------------- ##
  stopifnot(is.data.frame(coordinates))
  if (!all(c("chr", "pos") %in% names(coordinates))) {
    stop("'coordinates' must contain 'chr' and 'pos' columns")
  }
  if (!"id" %in% names(coordinates)) coordinates$id <- seq_len(nrow(coordinates))

  ## ---------------- GTF cache ---------------- ##
  gtf_result <- load_cached_gtf(gtffile, force = !cache_gtf)
  if (include_exon_info && (is.null(gtf_result$canonical_exons_by_gene) ||
                            length(gtf_result$canonical_exons_by_gene) == 0)) {
    message("[Cache rebuild] Canonical exon information missing; rebuilding cache…")
    gtf_result <- load_cached_gtf(gtffile, force = TRUE)
  }
  gtf_data  <- gtf_result$gtf_data
  genes     <- gtf_result$genes
  canon_by_gene <- gtf_result$canonical_exons_by_gene

  # NEW: Get canonical transcript IDs for each gene
  canon_tx_by_gene <- gtf_result$canonical_tx_by_gene

  ## ---------------- Chrom mapping ---------------- ##
  chr_mapping <- setNames(rep(NA_character_, length(unique(coordinates$chr))), unique(coordinates$chr))
  for (chr in names(chr_mapping)) {
    if (chr %in% seqlevels(genes))                     chr_mapping[chr] <- chr
    else if (sub("^chr", "", chr) %in% seqlevels(genes)) chr_mapping[chr] <- sub("^chr", "", chr)
    else if (paste0("chr", chr) %in% seqlevels(genes))    chr_mapping[chr] <- paste0("chr", chr)
    else warning("Chromosome ", chr, " absent from GTF; skipping positions on this contig.")
  }

  genes_by_chr <- lapply(seqlevels(genes), function(ch) genes[seqnames(genes) == ch])
  names(genes_by_chr) <- seqlevels(genes)

  ## ---------------- Results scaffold ---------------- ##
  n <- nrow(coordinates)
  results <- data.frame(
    id = coordinates$id,
    chr = coordinates$chr,
    position = coordinates$pos,
    genome = rep(genome, n),
    location_type = NA_character_,
    feature = NA_character_,
    gene_symbol = NA_character_,
    gene_strand = NA_character_,
    gene_type = NA_character_,
    upstream_gene = NA_character_,
    upstream_distance = NA_real_,
    downstream_gene = NA_character_,
    downstream_distance = NA_real_,
    within_exon = NA_integer_,
    fiveprime_exon = NA_integer_,
    threeprime_exon = NA_integer_,
    stringsAsFactors = FALSE
  )

  ## ---------------- Main loop ---------------- ##
  message("Annotating ", n, " positions…")
  for (i in seq_len(n)) {
    if (i %% 1000 == 0 || i == n) message("  processed ", i, "/", n)
    try({
      chr_raw <- coordinates$chr[i]
      pos_raw <- coordinates$pos[i]
      chr_use <- chr_mapping[[chr_raw]]
      if (is.na(chr_use)) next  # contig not in GTF

      query <- GRanges(chr_use, IRanges(pos_raw, pos_raw))
      chr_genes <- genes_by_chr[[chr_use]]
      if (!length(chr_genes)) next

      ov_gene <- findOverlaps(query, chr_genes, select = "first")
      if (!is.na(ov_gene)) {
        ## ---------------- Genic ---------------- ##
        results$location_type[i] <- "genic"
        g <- chr_genes[ov_gene]

        # FIXED: Properly access the symbol from mcols
        if ("symbol" %in% names(mcols(g))) {
          results$gene_symbol[i] <- mcols(g)$symbol
        } else if ("gene_name" %in% names(mcols(g))) {
          results$gene_symbol[i] <- mcols(g)$gene_name
        } else {
          results$gene_symbol[i] <- mcols(g)$gene_id
        }

        results$gene_strand[i] <- as.character(strand(g))
        results$gene_type[i] <- if ("gene_type" %in% names(mcols(g))) mcols(g)$gene_type
        else if ("gene_biotype" %in% names(mcols(g))) mcols(g)$gene_biotype
        else NA_character_

        # Process feature information
        gene_id <- mcols(g)$gene_id

        # NEW: Get canonical transcript ID for this gene
        canon_tx_id <- NULL
        if (!is.null(canon_tx_by_gene) && gene_id %in% names(canon_tx_by_gene)) {
          canon_tx_id <- canon_tx_by_gene[[gene_id]]
        }

        # MODIFIED: Only get exons for the canonical transcript
        exons_g <- NULL
        if (!is.null(canon_tx_id)) {
          exons_g <- gtf_data[gtf_data$type == "exon" &
                                gtf_data$gene_id == gene_id &
                                gtf_data$transcript_id == canon_tx_id]
        } else {
          # Fall back to all exons if no canonical transcript is defined
          exons_g <- gtf_data[gtf_data$type == "exon" & gtf_data$gene_id == gene_id]
        }

        utr5_g <- gtf_data[gtf_data$type == "five_prime_utr" & gtf_data$gene_id == gene_id]
        if (!is.null(canon_tx_id)) {
          utr5_g <- utr5_g[utr5_g$transcript_id == canon_tx_id]
        }

        utr3_g <- gtf_data[gtf_data$type == "three_prime_utr" & gtf_data$gene_id == gene_id]
        if (!is.null(canon_tx_id)) {
          utr3_g <- utr3_g[utr3_g$transcript_id == canon_tx_id]
        }

        feat <- NULL
        if (length(utr5_g) && !is.na(findOverlaps(query, utr5_g, select = "first"))) {
          feat <- "5' UTR"
        } else if (length(utr3_g) && !is.na(findOverlaps(query, utr3_g, select = "first"))) {
          feat <- "3' UTR"
        } else if (length(exons_g)) {
          ov_ex <- findOverlaps(query, exons_g, select = "first")
          if (!is.na(ov_ex)) {
            ex_hit <- exons_g[ov_ex]
            if ("exon_number" %in% names(mcols(ex_hit))) {
              # NEW: Use canonical exon numbering for feature
              ex_num <- as.integer(mcols(ex_hit)$exon_number)
              feat <- paste0("exon ", ex_num)
              results$within_exon[i] <- ex_num  # Set within_exon to match feature
            } else {
              feat <- "exon"
            }
          }
        }
        results$feature[i] <- if (!is.null(feat)) feat else "intron"

        ## ----- Canonical‑exon extras ----- ##
        if (include_exon_info) {
          exons_can <- canon_by_gene[[gene_id]]
          if (!is.null(exons_can) && length(exons_can)) {
            ov_can <- findOverlaps(query, exons_can, select = "first")
            if (!is.na(ov_can)) {
              # Should already be set above if we're using consistent exon numbering
              if (is.na(results$within_exon[i])) {
                results$within_exon[i] <- as.integer(mcols(exons_can[ov_can])$exon_number)
              }
            } else {
              ## intronic: nearest exons
              ends_can   <- end(exons_can)
              starts_can <- start(exons_can)
              up_idx  <- which(ends_can < pos_raw)
              dn_idx  <- which(starts_can > pos_raw)
              if (length(up_idx)) {
                results$fiveprime_exon[i] <- as.integer(mcols(exons_can[up_idx[which.max(ends_can[up_idx])]])$exon_number)
              }
              if (length(dn_idx)) {
                results$threeprime_exon[i] <- as.integer(mcols(exons_can[dn_idx[which.min(starts_can[dn_idx])]])$exon_number)
              }
            }
          }
        }
      } else {
        ## ---------------- Intergenic ---------------- ##
        results$location_type[i] <- "intergenic"
        up <- chr_genes[start(chr_genes) > pos_raw]
        dn <- chr_genes[end(chr_genes) < pos_raw]

        if (length(up)) {
          d_up <- start(up) - pos_raw
          idx <- which.min(d_up)
          # FIXED: Properly access the symbol from mcols
          if ("symbol" %in% names(mcols(up))) {
            results$upstream_gene[i] <- mcols(up)$symbol[idx]
          } else {
            results$upstream_gene[i] <- mcols(up)$gene_id[idx]
          }
          results$upstream_distance[i] <- d_up[idx]
        }

        if (length(dn)) {
          d_dn <- pos_raw - end(dn)
          idx <- which.min(d_dn)
          # FIXED: Properly access the symbol from mcols
          if ("symbol" %in% names(mcols(dn))) {
            results$downstream_gene[i] <- mcols(dn)$symbol[idx]
          } else {
            results$downstream_gene[i] <- mcols(dn)$gene_id[idx]
          }
          results$downstream_distance[i] <- d_dn[idx]
        }
      }
    }, silent = TRUE)
  }
  message("Annotation complete.")
  results
}


#' Preload GTF data with proper symbol handling and canonical transcript tracking
#' @export
#' @keywords internal
preload_gtf <- function(gtffile) {
  suppressPackageStartupMessages({ library(digest); library(rtracklayer); library(data.table); library(GenomicRanges) })
  message("[GTF preload] ", basename(gtffile))
  if (!exists("gtf_env", .GlobalEnv)) assign("gtf_env", new.env(), .GlobalEnv)
  env <- get("gtf_env", .GlobalEnv); key <- paste0("gtf_", digest(gtffile))

  ## Import GTF data
  message("[Importing GTF] ", basename(gtffile))
  gtf <- import(gtffile, format="gtf")
  message("[Processing Gene Level Data] ", basename(gtffile))
  keep <- gtf$type %in% c("gene","transcript","exon","five_prime_utr","three_prime_utr")
  gtf <- gtf[keep]

  ## Process gene data with proper symbol handling
  genes <- gtf[gtf$type == "gene"]

  # Convert to data frame to easily check column existence and properly handle NAs
  genes_df <- as.data.frame(mcols(genes))

  # Add symbol column based on available gene identifiers
  if ("gene_name" %in% colnames(genes_df)) {
    mcols(genes)$symbol <- ifelse(!is.na(genes_df$gene_name), genes_df$gene_name, genes_df$gene_id)
  } else {
    # Fall back to gene_id if gene_name doesn't exist
    mcols(genes)$symbol <- genes_df$gene_id
  }

  ## Process transcript data
  message("[Processing Transcript Level Data] ", basename(gtffile))
  transcripts <- gtf[gtf$type == "transcript" & !is.na(gtf$gene_id)]

  # NEW: Store canonical transcript IDs by gene
  canon_tx_by_gene <- list()

  ## Process canonical transcript and exon data
  message("[Processing Exon Level Data] ", basename(gtffile))
  canon_by_gene <- list()

  if (length(transcripts) > 0 && "transcript_id" %in% names(mcols(transcripts))) {
    # Identify canonical (longest) transcript for each gene
    tx_dt <- as.data.table(as.data.frame(transcripts)[, c("gene_id","transcript_id","start","end")])
    tx_dt[, len := end - start + 1]

    # Get canonical transcript ID for each gene
    canon_tx_dt <- tx_dt[, .SD[which.max(len)], by=gene_id]

    # Create gene_id to canonical transcript_id mapping
    for (i in 1:nrow(canon_tx_dt)) {
      gene_id <- canon_tx_dt$gene_id[i]
      tx_id <- canon_tx_dt$transcript_id[i]
      canon_tx_by_gene[[gene_id]] <- tx_id
    }

    # Get canonical exons
    exons_can <- gtf[gtf$type == "exon" & gtf$transcript_id %in% canon_tx_dt$transcript_id]
    if (length(exons_can) > 0) {
      canon_by_gene <- split(exons_can, exons_can$gene_id, drop=TRUE)
    }
  }

  ## Cache results
  assign(key, gtf, envir=env)
  assign(paste0(key,"_genes"), genes, envir=env)
  assign(paste0(key,"_canonical_exons"), canon_by_gene, envir=env)
  assign(paste0(key,"_canonical_tx"), canon_tx_by_gene, envir=env)

  message("[GTF Preload Complete] ", basename(gtffile))
  invisible(TRUE)
}

#' Load cached GTF data or rebuild if needed
#' @export
#' @keywords internal
load_cached_gtf <- function(gtffile, force = FALSE) {
  suppressPackageStartupMessages({ library(digest) })
  if (!exists("gtf_env", .GlobalEnv)) assign("gtf_env", new.env(), .GlobalEnv)
  env <- get("gtf_env", .GlobalEnv); key <- paste0("gtf_", digest(gtffile))
  if (force || !exists(key, envir=env)) preload_gtf(gtffile)

  # Get canonical transcript IDs for each gene (new)
  canon_tx_key <- paste0(key,"_canonical_tx")
  canon_tx_by_gene <- if (exists(canon_tx_key, envir=env)) get(canon_tx_key, envir=env) else list()

  list(
    gtf_data = get(key, envir=env),
    genes = get(paste0(key,"_genes"), envir=env),
    canonical_exons_by_gene = get(paste0(key,"_canonical_exons"), envir=env),
    canonical_tx_by_gene = canon_tx_by_gene
  )
}



#' Convenience function to annotate a single position
#' @export
#' @keywords internal
annotate_genomic_position <- function(chr, pos, genome, gtffile,
                                      tss_upstream=2000, tss_downstream=200,
                                      cache_gtf=TRUE, include_exon_info=FALSE) {
  annotate_genomic_coordinates(
    data.frame(chr=chr, pos=pos),
    genome, gtffile, tss_upstream, tss_downstream, cache_gtf, include_exon_info
  )[1,]
}

