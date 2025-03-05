roxygen2::roxygenize()

data_S1 <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL1_S1/20250302_232026_humanwgs_singleton/outputs.json")
data_S18 <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL2_S18/20250302_232206_humanwgs_singleton/outputs.json")
data_S42 <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL3_S42/20250303_062244_humanwgs_singleton/outputs.json")
data_SAdd <- read_json_file("/Volumes/furlan_s/sfurlan/250302_leuklong/pbWGS/LL4_SAdd/20250302_231552_humanwgs_singleton/outputs.json")

file <- "/Users/sfurlan/Library/CloudStorage/OneDrive-SharedLibraries-FredHutchinsonCancerResearchCenter/Furlan_Lab - General/experiments/leuklong/250302_SCRI/res/LL4_SAdd.small_variants.ann.txt"

plot_snpEff(file)

anns <- vroom(file, guess_max = 1000000)
flt3 <- anns[(anns$CHROM %in% "chr13" & anns$POS > 28003274 & anns$POS < 28100592),]
wt1 <- anns[(anns$CHROM %in% "chr11" & anns$POS > 32396330 & anns$POS <  32396400),]

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
