# Load required packages
library(GenomicRanges)
library(rtracklayer)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Retrieve gene annotations from TxDb (GRCh38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genesGR <- genes(txdb)

# Annotate genes with symbols using org.Hs.eg.db
gene_ids <- names(genesGR)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
mcols(genesGR)$symbol <- gene_symbols

# Define a helper function to extend GRanges (similar to ArchR's extendGR)
extendGR <- function(gr, upstream, downstream) {
  new_gr <- gr
  # For genes on the positive strand (or no strand), extend start by upstream and end by downstream.
  pos_idx <- which(as.character(strand(gr)) %in% c("+", "*"))
  # For genes on the negative strand, reverse: extend start by downstream and end by upstream.
  neg_idx <- which(as.character(strand(gr)) == "-")

  if(length(pos_idx) > 0){
    start(new_gr)[pos_idx] <- pmax(1, start(gr)[pos_idx] - upstream)
    end(new_gr)[pos_idx] <- end(gr)[pos_idx] + downstream
  }
  if(length(neg_idx) > 0){
    start(new_gr)[neg_idx] <- pmax(1, start(gr)[neg_idx] - downstream)
    end(new_gr)[neg_idx] <- end(gr)[neg_idx] + upstream
  }
  new_gr
}

# Default parameters from the addGeneScoreMatrix function:
geneUpstream <- 5000  # extend gene body upstream
geneDownstream <- 0   # no extension downstream
useTSS <- FALSE       # use full gene body model (not TSS only)
extendTSS <- TRUE     # if useTSS were TRUE, extend the TSS

# Extend gene regions according to the defaults:
if(!useTSS){
  # When useTSS is FALSE, extend the full gene body.
  extendedGenes <- extendGR(genesGR, upstream = geneUpstream, downstream = geneDownstream)
} else {
  # If using TSS, first shrink the gene to a 1bp region at the TSS (for + strand, start; for - strand, end)
  extendedGenes <- genesGR
  pos_idx <- which(as.character(strand(extendedGenes)) %in% c("+", "*"))
  neg_idx <- which(as.character(strand(extendedGenes)) == "-")
  if(length(pos_idx) > 0){
    end(extendedGenes)[pos_idx] <- start(extendedGenes)[pos_idx]
  }
  if(length(neg_idx) > 0){
    start(extendedGenes)[neg_idx] <- end(extendedGenes)[neg_idx]
  }
  # Then extend these TSS regions if desired.
  if(extendTSS){
    extendedGenes <- extendGR(extendedGenes, upstream = geneUpstream, downstream = geneDownstream)
  }
}

# Write the extended gene regions to a BED file
export(extendedGenes, file.path("/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac", "GRCh38_extendedGenes.bed"), format = "BED")

