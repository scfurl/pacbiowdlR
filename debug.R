roxygen2::roxygenize()
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(rtracklayer)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genesGR <- genes(txdb)
# # Annotate genes with symbols using org.Hs.eg.db
gene_ids <- names(genesGR)
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")
mcols(genesGR)$name <- gene_symbols
mcols(genesGR)$symbol <- gene_symbols
#debug(generateExtendedGeneIntervals)
genesExtended <- generateExtendedGeneIntervals(genesGR, cores = parallel::detectCores())
extendedDF <- data.frame(
  seqnames = seqnames(genesExtended),
  start = start(genesExtended) - 1,  # BED format is 0-based
  end = end(genesExtended),
  name = mcols(genesExtended)$symbol,
  weight = mcols(genesExtended)$geneWeight,
  strand = strand(genesExtended),
  tss = mcols(genesExtended)$tss,
  geneStart = mcols(genesExtended)$geneStart,
  geneEnd = mcols(genesExtended)$geneEnd
)

outputFile <- file.path("/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac", "GRCh38_extended.bed")
write.table(
  extendedDF,
  file = outputFile,
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)


bed_df <- data.frame(
  chrom = seqnames(genesGR),
  start = start(genesGR) - 1,  # BED format is 0-based for start positions
  end = end(genesGR),
  name = mcols(genesGR)$name,
  strand = as.character(strand(genesGR))
)

write.table(bed_df, file.path("/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac", "GRCh38.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

















bed_df <- data.frame(
  chrom = seqnames(genesGR),
  start = start(genesGR) - 1,  # BED format is 0-based for start positions
  end = end(genesGR),
  name = mcols(genesGR)$name,
  strand = as.character(strand(genesGR))
)

write.table(bed_df, file.path("/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac", "GRCh38.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
#xsdebug(createGeneRegions)
windows <- createGeneRegions(genes = genesGR, includeChr = paste0("chr", c(1:22, "X")))

windows_df <- data.frame(
  chrom = seqnames(windows),
  start = start(windows) - 1,  # BED format is 0-based for start positions
  end = end(windows),
  name = mcols(windows)$name,
  strand = as.character(strand(windows))
)

write.table(windows_df, file.path("/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerCenter/Furlan_Lab - General/experiments/leuklong/integration/fire_atac", "GRCh38_extended.bed"), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

# go and run bedtools map on cluster to obtain FIRE score over these extended gene regions

ROOT_DIR2 <- "/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerResearchCenter/Furlan_Lab - General/datasets/Healthy_BM_greenleaf"
seur <- readRDS(file.path(ROOT_DIR2, "230329_rnaAugmented_seurat.RDS"))

library(Seurat)
library(flscuts)
library(monocle3)
library(patchwork)
DimPlot(seur, group.by = "SFClassification", cols = seur@misc$colors)

cds <- seurat_to_monocle3(seur)



data_S1 <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL1_S1/20250302_232026_humanwgs_singleton/outputs.json")
data_S18 <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL2_S18/20250302_232206_humanwgs_singleton/outputs.json")
data_S42 <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL3_S42/20250303_062244_humanwgs_singleton/outputs.json")
data_SAdd <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL4_SAdd/20250302_231552_humanwgs_singleton/outputs.json")

file <- "/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerResearchCenter/Furlan_Lab - General/experiments/leuklong/250302_SCRI/res/LL4_SAdd.small_variants.ann.txt"

DATA_DIR <- "/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerResearchCenter/Furlan_Lab - General/experiments/leuklong/integration"
files <- list.files(file.path(DATA_DIR, "fire_atac"))
file <- files[grepl("GA", files)]
files <- file.path(DATA_DIR, "fire_atac", file)
samps <- c( "EP300-ZNF384_ALL", "ETV6-RUNX1_ALL", "APML","NUP-NSD1_AML")

make_matrix_from_fs_gas <- function(files, samps=NULL){
 if(is.null(samps)){
   samps <- basename(files)
 }
 if(length(samps)!=length(files)) stop("length of files vector and samps should be the same")
 bed_data <- lapply(1:length(files), function(x) {
   dat <- data.table::fread(files[x], header = F)
   dat$samp <- samps[x]
   colnames(dat)<-c("seqnames", "start", "end", "gene", "strand", "score", "samp")
   dat$score[dat$score=="."] <- 0
   dat$score[dat$score < 0] <- 0
   gas_vector <- dat$score
   gas_vector <- as.numeric(gas_vector)
   names(gas_vector) <- dat$gene
   gas_vector
 })
 mat <- Matrix(do.call(cbind, bed_data), sparse = T)
 colnames(mat)<-samps
 mat
}

#debug(make_matrix_from_fs_gas)
mat <- make_matrix_from_fs_gas(files, samps)




fs <- SummarizedExperiment(assays = as.matrix(mat), colData = DataFrame(name = colnames(mat), type = colnames(mat), row.names = colnames(mat)))

PCs <- 25 #Number of PCs for clustering
n_top <- c(3000,2500,2500) #number of features
rownames(seur[["umap"]]@cell.embeddings)


rna<-iterative_LSI(cds, num_dim = PCs, num_features = n_top, resolution = c(2e-4, 9e-4, 1e-3), verbose = T)
rna <- reduce_dimension(rna, reduction_method = "UMAP", preprocess_method = "LSI", verbose = T, umap.min_dist = 0.45, umap.n_neighbors = 55,umap.metric = "euclidean")

p1 <- plot_cells(rna, color_cells_by = "SFClassification")+scale_color_manual(values = seur@misc$colors)

projection <- project_data(projector = rna, projectee = fs, reduced_dim = "LSI", embedding = "UMAP", features = "annotation-based", projectee_label_col = "type", make_pseudo_single_cells = T, n = 1000, force = T)


p<-as.data.frame(do.call(rbind, projection[2:1]))
cols2 <- pals::glasbey(21)[3:(3+dim(fs)[2])]
names(cols2)<-colData(fs)$name
cols2 <- c(cols2, "single_cell_reference"="lightgray")

p2<-ggplot(p, aes(x=UMAP1, y=UMAP2, color=projectee_labels))+geom_point(size=0.2, alpha=0.5)+ guides(colour = guide_legend(override.aes = list(size=5)))+monocle3:::monocle_theme_opts()+scale_color_manual(values=cols2 )
p1+p2


fs




muts <- as.data.frame(rbind(wt1[1,], flt3[51,])) %>% as.data.frame()
write.csv(muts, "~/Desktop/muts.csv")

annotated_vcf_table <- "/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerResearchCenter/Furlan_Lab - General/experiments/leuklong/250302_SCRI/res/LL4_SAdd.small_variants.clinvar.txt"

ann <- vroom(annotated_vcf_table, guess_max = 1000000) |>
  filter(!is.na(ALLELEID))
patho <- filter(ann, grepl("Pathogenic",CLNSIG) | grepl("Likely_pathogenic",CLNSIG)) |> separate("ID", c("ID","unknown"), sep = ";")
#unique(patho$CLNSIG)
#nrow(patho)
as.data.frame(patho) |> gt()

flt3 <- ann[grep("^FLT3", ann$GENEINFO),]
as.data.frame(flt3) |> gt()
wt1 <- ann[grep("^WT1", ann$GENEINFO),]
as.data.frame(wt1) |> gt()




plot_circos_from_vcf(data_S42$humanwgs_singleton.tertiary_sv_filtered_vcf, thresh = 15, return_data = T)
plot_circos_from_vcf(data_SAdd$humanwgs_singleton.tertiary_sv_filtered_vcf, highlight = 12, return_data = F, title = "NUP98::NSD1 AML")
plot_circos_from_vcf(data_S1$humanwgs_singleton.tertiary_sv_filtered_vcf, highlight = 3, return_data = F, title = "EP300::ZNF384 ALL")
plot_circos_from_vcf(data_S18$humanwgs_singleton.tertiary_sv_filtered_vcf, highlight = 3, return_data = F, thresh = 15, title = "ETV6-RUNX1 ALL")
plot_circos_from_vcf(data_S42$humanwgs_singleton.tertiary_sv_filtered_vcf, highlight = 5, return_data = F, thresh = 15, title = "PML::RARA APML")






#BiocManager::install("org.Hs.eg.db", force=T)
library(pacbiowdlR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

  # Assign colors to chromosomes
#chrom_palette <- scales::hue_pal()(24)  # 23 distinct colors

chrom_palette <- c(
  "#FF0000", "#FF9900", "#FFCC00", "#00FF00", "#6699FF", "#CC33FF",
  s"#999912", "#999999", "#FF00CC", "#CC0000", "#FFCCCC", "#FFFF00",
  "#CCFF00", "#358000", "#0000CC", "#99CCFF", "#00FFFF", "#ECFFFF",
  "#9900CC", "#CC99FF", "#996600", "#666600", "#666666", "#CCCCCC",
  "#79CC3B", "#E0EC3B", "#CCC99F"
)
chr_colors <- setNames(chrom_palette, paste0("chr", c(1:22, "X")))


generateCNAPlotDiscoverGenes2(
  depth_bigwig_file  =   data_S1$humanwgs_singleton.cnv_depth_bw,
  variant_file       =   data_S1$humanwgs_singleton.cnv_vcf,
  txdb               = TxDb.Hsapiens.UCSC.hg38.knownGene,
  downsample         = 0.1,
  samplename         = "MySample", gene_delta_threshold = 500,
  return_data        = F, colors = chrom_palette
)

generateCNAPlotDiscoverGenes(
  depth_bigwig_file  =   data_S1$humanwgs_singleton.cnv_depth_bw,
  variant_file       =   data_S1$humanwgs_singleton.cnv_vcf,
  txdb               = TxDb.Hsapiens.UCSC.hg38.knownGene,
  downsample         = 0.8,
  samplename         = "MySample", highlight_genes = c("CDKN2A", "CDKN2B", "p14ARF", "MLLT3",  "MTAP"),
  return_data        = F, colors = chrom_palette, chr_filter = "chr9"
)



generateCNAPlotDiscoverGenes2(
  depth_bigwig_file  =   data_S18$humanwgs_singleton.cnv_depth_bw,
  variant_file       =   data_S18$humanwgs_singleton.cnv_vcf,
  txdb               = TxDb.Hsapiens.UCSC.hg38.knownGene,
  downsample         = 0.01,
  samplename         = "MySample", gene_delta_threshold = 500,
  return_data        = F, colors = chrom_palette
)

generateCNAPlotDiscoverGenes(
  depth_bigwig_file  =   data_S18$humanwgs_singleton.cnv_depth_bw,
  variant_file       =   data_S18$humanwgs_singleton.cnv_vcf,
  txdb               = TxDb.Hsapiens.UCSC.hg38.knownGene,
  downsample         = 0.8,
  samplename         = "MySample", highlight_genes = "ETV6",
  return_data        = F, colors = chrom_palette, chr_filter = "chr12"
)

generateCNAPlotDiscoverGenes(
  depth_bigwig_file  =   data_S18$humanwgs_singleton.cnv_depth_bw,
  variant_file       =   data_S18$humanwgs_singleton.cnv_vcf,
  txdb               = TxDb.Hsapiens.UCSC.hg38.knownGene,
  downsample         = 0.8,
  samplename         = "MySample", highlight_genes = "RUNX1",
  return_data        = F, colors = chrom_palette, chr_filter = "chr21"
)


generateCNAPlotDiscoverGenes(
  depth_bigwig_file  =   data_S42$humanwgs_singleton.cnv_depth_bw,
  variant_file       =   data_S42$humanwgs_singleton.cnv_vcf,
  txdb               = TxDb.Hsapiens.UCSC.hg38.knownGene,
  downsample         = 0.01,
  samplename         = "MySample", gene_delta_threshold = 500,
  return_data        = F, colors = chrom_palette
)

generateCNAPlotDiscoverGenes(
  depth_bigwig_file  =   data_SAdd$humanwgs_singleton.cnv_depth_bw,
  variant_file       =   data_SAdd$humanwgs_singleton.cnv_vcf,
  txdb               = TxDb.Hsapiens.UCSC.hg38.knownGene,
  downsample         = 0.01,
  samplename         = "MySample", gene_delta_threshold = 500,
  return_data        = F, colors = chrom_palette
)

d <- jsonlite::read_json(data_SAdd$humanwgs_singleton.pharmcat_report_json)
names(d$genes[[1]])

readLines(data_SAdd$humanwgs_singleton.pharmcat_report_html)

d$results[[2]]

sv <- read.table(data_SAdd$humanwgs_singleton.tertiary_sv_filtered_tsv, sep = "\t")
head(sv)
sv <- sv[sv$V6=="BND",]

tabb<-pharmcat_json(data_S1$humanwgs_singleton.pharmcat_report_json)
gt_tbl <- gt::gt(tabb)

write.csv(tabb, "~/Desktop/pharmcat.csv")
# Show the gt Table
gt_tbl



BiocManager::install("chromVAR")
devtools::install_github("GreenleafLab/ArchR", ref="master", repos = BiocManager::repositories())
library(ArchR)
ArchR::installExtraPackages()
