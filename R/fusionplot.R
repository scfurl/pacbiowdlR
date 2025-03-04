# Load required libraries
library(circlize)
library(VariantAnnotation)
library(dplyr)

# Function to read gzipped VCF and generate a Circos plot
plot_circos_from_vcf <- function(vcf_file,title = "", thresh = 20, allowed = paste0("chr", c(1:22)), highlight = NULL, return_data = T){

  # Check if file exists
  if (!file.exists(vcf_file)) {
    stop(paste("Error: File not found:", vcf_file))
  }

  # Read VCF file
  vcf <- readVcf(vcf_file, genome = "GRCh38")
  vcf <- vcf[vcf@info$SVTYPE=="BND",]
  vcf <- vcf[vcf@fixed$FILTER=="PASS"]
  vcf <- vcf[vcf@info$NotFullySpanned=="FALSE"]
  vcf <- vcf[grepl("protein_coding", sapply(vcf@info$BCSQ, function(bs){ifelse(length(bs)==0, NA, bs[1])})),]
  vcf <- vcf[vcf@assays@data$DP[,1]>thresh,]

  dat <- gsub("pbsv.BND.", "", as.character(info(vcf)$MATEID))
  parts <- strsplit(dat, "-")
  sv_data <- do.call(rbind, parts)
  sv_data <-data.frame(chr1=strsplit(sv_data[,1], ":") %>% sapply("[[", 1),
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

  if(return_data){
    return(sv_data)
  } else {
    highlight <- sv_data[highlight, ]
    bed3<-highlight[,c(1,3)]
    colnames(bed3) <- c("chr", "start")
    bed3$start = as.numeric(bed3$start )
    bed3$end <- bed3$start
    bed4 = as.data.frame(highlight[,c(2,4)])
    colnames(bed4) <- c("chr", "start")
    bed4$start = as.numeric(bed4$start )
    bed4$end <- bed4$start

    bed1 = sv_data[,c(1,3)]
    colnames(bed1) <- c("chr", "start")
    bed1$start = as.numeric(bed1$start )
    bed1$end <- bed1$start
    bed2 = as.data.frame(sv_data[,c(2,4)])
    colnames(bed2) <- c("chr", "start")
    bed2$start = as.numeric(bed2$start )
    bed2$end <- bed2$start


    circos.initializeWithIdeogram()
    circos.genomicLink(bed1, bed2, col = rand_color(nrow(bed1), transparency = 0.5),
                       border = NA)
    circos.genomicLink(bed3, bed4, col = "black",
                       lwd = 10)
    circos.clear()
    title(title)
  }


}

