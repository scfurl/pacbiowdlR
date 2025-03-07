# # Load required packages
# library(GenomicRanges)
# library(rtracklayer)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
# library(org.Hs.eg.db)
#
# # Retrieve gene annotations from TxDb (GRCh38)
# txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
# genesGR <- genes(txdb)
#
# # Annotate genes with symbols using org.Hs.eg.db
# gene_ids <- names(genesGR)
# gene_symbols <- mapIds(org.Hs.eg.db,
#                        keys = gene_ids,
#                        column = "SYMBOL",
#                        keytype = "ENTREZID",
#                        multiVals = "first")
# mcols(genesGR)$symbol <- gene_symbols
#
# # Define a helper function to extend GRanges (similar to ArchR's extendGR)
# extendGR <- function(gr, upstream, downstream) {
#   new_gr <- gr
#   # For genes on the positive strand (or no strand), extend start by upstream and end by downstream.
#   pos_idx <- which(as.character(strand(gr)) %in% c("+", "*"))
#   # For genes on the negative strand, reverse: extend start by downstream and end by upstream.
#   neg_idx <- which(as.character(strand(gr)) == "-")
#
#   if(length(pos_idx) > 0){
#     start(new_gr)[pos_idx] <- pmax(1, start(gr)[pos_idx] - upstream)
#     end(new_gr)[pos_idx] <- end(gr)[pos_idx] + downstream
#   }
#   if(length(neg_idx) > 0){
#     start(new_gr)[neg_idx] <- pmax(1, start(gr)[neg_idx] - downstream)
#     end(new_gr)[neg_idx] <- end(gr)[neg_idx] + upstream
#   }
#   new_gr
# }
#
# # Default parameters from the addGeneScoreMatrix function:
# geneUpstream <- 5000  # extend gene body upstream
# geneDownstream <- 0   # no extension downstream
# useTSS <- FALSE       # use full gene body model (not TSS only)
# extendTSS <- TRUE     # if useTSS were TRUE, extend the TSS
#
# # Extend gene regions according to the defaults:
# if(!useTSS){
#   # When useTSS is FALSE, extend the full gene body.
#   extendedGenes <- extendGR(genesGR, upstream = geneUpstream, downstream = geneDownstream)
# } else {
#   # If using TSS, first shrink the gene to a 1bp region at the TSS (for + strand, start; for - strand, end)
#   extendedGenes <- genesGR
#   pos_idx <- which(as.character(strand(extendedGenes)) %in% c("+", "*"))
#   neg_idx <- which(as.character(strand(extendedGenes)) == "-")
#   if(length(pos_idx) > 0){
#     end(extendedGenes)[pos_idx] <- start(extendedGenes)[pos_idx]
#   }
#   if(length(neg_idx) > 0){
#     start(extendedGenes)[neg_idx] <- end(extendedGenes)[neg_idx]
#   }
#   # Then extend these TSS regions if desired.
#   if(extendTSS){
#     extendedGenes <- extendGR(extendedGenes, upstream = geneUpstream, downstream = geneDownstream)
#   }
# }
#
#
# # Write the extended gene regions to a BED file
# export(extendedGenes, file.path("/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac", "GRCh38_extendedGenes.bed"), format = "BED")
#
# library(GenomicRanges)
# library(rtracklayer)
#
# # Define a default gene model function (adjust as needed)
# defaultGeneModel <- function(x) {
#   exp(-abs(x)/5000) + exp(-1)
# }
#
# library(doParallel)
# library(foreach)
#
#
#
# calculateFiberSeqFIREScore <- function(fire_file, defaultGeneModel, cores = 8) {
#   # If fire is a character, assume it's a file path to a BED file and import it.
#   if (is.character(fire_file)) {
#     message("Importing FIRE BED file...")
#     fire <- data.table::fread(fire_file)
#     colnames(fire) <- c("seqnames", "start", "end", "gene_id", "something", "strand", "score")
#     fire <- GRanges(fire)
#     fire$score[fire$score=="."]<-"0"
#     fire$score <- as.numeric(fire$score)
#     #fire$symbol <- as.character(extendedGenes$symbol[match(fire$gene_id, extendedGenes$gene_id)])
#   }
#
#   # Ensure the 'fire' GRanges has a metadata column named "score"
#   if (!"score" %in% colnames(mcols(fire))) {
#     stop("The FIRE GRanges must have a 'score' column in its metadata.")
#   }
#
#
#
#   cl <- makeCluster(cores)
#   registerDoParallel(cl)
#
#   # Parallel loop using foreach
#   geneScores <- foreach(g = seq_along(fire),
#                         .combine = c,
#                         .packages = "GenomicRanges") %dopar% {
#                           peakIdx <- subjectHits(hits)[queryHits(hits) == g]
#                           if (length(peakIdx) == 0) {
#                             0
#                           } else {
#                             # Use the center of the gene as reference (using original gene width)
#                             geneCenter <- start(extendedGenes[g]) + (width(genes[g]) / 2)
#                             # Calculate peak centers
#                             peakCenters <- start(fire[peakIdx]) + (width(fire[peakIdx]) / 2)
#                             distances <- peakCenters - geneCenter
#
#                             # Compute weights using the gene model function
#                             weights <- gene_model(distances)
#                             # Multiply weights by FIRE score (assumed to be in mcols(fire)$score)
#                             scores <- mcols(fire[peakIdx])$score * weights
#                             sum(scores, na.rm = TRUE)
#                           }
#                         }
#
#   # Stop the cluster
#   stopCluster(cl)
#
#   # Check the result
#   geneScores
#   return(geneScores)
# }
#
# # Example usage:
# # Assume 'genesGR' is your GRanges object for genes with a "symbol" column.
# # And 'fire_file' is the path to your FIRE BED file.
# # genesGR <- ... (load or create your gene GRanges)
# # fire_file <- "path/to/your/FIRE.bed"
# fire_file <- "/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerResearchCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac/LL4_SAdd_GA_sum_scores.bed"
# #fire <- "/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerResearchCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac/LL4_SAdd_GA_sum_scores.bed"
# geneFIREscores <- calculateFiberSeqFIREScore(fire_file)
# head(geneFIREscores)
#
#
#
#
