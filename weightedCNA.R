library(Rsamtools)      # for FASTA access
library(Biostrings)     # for sequence handling and letterFrequency()
library(GenomicRanges)  # for GRanges and Views
library(progress)       # for a simple progress bar

compute_interval_metrics_vectorized_progress <- function(fasta_file, intervals, repeat_gr, chunk_size = 100) {
  # Open the FASTA file
  fasta <- FaFile(fasta_file)
  open(fasta)

  n_intervals <- length(intervals)
  # Create a progress bar
  pb <- progress_bar$new(
    format = "  Processing [:bar] :percent ETA: :eta",
    total = n_intervals, clear = FALSE, width= 60
  )

  # Vectorize GC calculation: get sequences for all intervals at once if possible.
  # But here, we process in chunks to update progress.
  # Split intervals into chunks
  chunks <- split(intervals, ceiling(seq_along(intervals) / chunk_size))
  result_list <- vector("list", length(chunks))

  # Pre-calculate the repeat coverage profile once
  rep_cov <- coverage(repeat_gr)

  processed <- 0
  for (i in seq_along(chunks)) {
    chunk <- chunks[[i]]

    # Get sequences for the current chunk (vectorized)
    sequences <- getSeq(fasta, chunk)
    # Calculate GC content using a vectorized letterFrequency call
    gc_counts <- rowSums(letterFrequency(sequences, letters = c("G", "C"), as.prob = FALSE))
    widths <- width(sequences)
    gc_content <- gc_counts / widths

    # Calculate repeat fraction for each interval in the chunk using Views
    rep_frac <- sapply(seq_along(chunk), function(j) {
      interval <- chunk[j]
      chr <- as.character(seqnames(interval))
      if (!chr %in% names(rep_cov)) return(0)
      # Create a Views object on the repeat coverage for the region.
      v <- Views(rep_cov[[chr]], start = start(interval), end = end(interval))
      # Count the number of bases with a nonzero repeat annotation.
      n_rep <- sum(as.numeric(unlist(v)) > 0)
      n_rep / width(interval)
    })

    # Construct a data frame for the chunk
    chunk_df <- data.frame(
      seqnames = as.character(seqnames(chunk)),
      start = start(chunk),
      end = end(chunk),
      width = widths,
      gc_content = gc_content,
      repeat_fraction = rep_frac,
      stringsAsFactors = FALSE
    )

    result_list[[i]] <- chunk_df
    processed <- processed + length(chunk)
    pb$update(processed / n_intervals)   # update progress bar
  }

  close(fasta)
  result <- do.call(rbind, result_list)
  return(result)
}





# (1) Load your intervals from your BigWig file.
coverage_data <- import.bw(data_4N$humanwgs_singleton.cnv_depth_bw)
coverage_data <- keepStandardChromosomes(coverage_data, pruning.mode = "tidy")
#coverage_data <- coverage_data[sample(1:length(coverage_data), 1419),]

# (2) Process your RepeatMasker file to get repeat annotations.
repeat_df <- read.table("/Users/sfurlan/refs/GRCh38/hg38.fa.out",
                        skip = 2,
                        fill = TRUE,
                        header = FALSE,
                        stringsAsFactors = FALSE)
colnames(repeat_df) <- c("SW_score", "perc_div", "perc_del", "perc_ins",
                         "query_sequence", "query_begin", "query_end", "query_left",
                         "strand", "repeat", "repeat_class_family",
                         "repeat_begin", "repeat_end", "repeat_left", "ID")
repeat_df$query_begin <- as.numeric(repeat_df$query_begin)
repeat_df$query_end <- as.numeric(repeat_df$query_end)
repeat_df$strand <- gsub("C", "*", repeat_df$strand)
repeat_gr <- GRanges(seqnames = repeat_df$query_sequence,
                     ranges = IRanges(start = repeat_df$query_begin,
                                      end = repeat_df$query_end),
                     strand = repeat_df$strand)

# (3) Set the path to your reference genome FASTA.
fasta_path <- "/Users/sfurlan/refs/GRCh38/GRCh38.p13.genome.fa"

# (4) For progress reporting, we use FutureParam instead of MulticoreParam.

# Set the future plan. For example, on Mac/Linux:
plan(multisession, workers = 16)

# Choose a progress handler.
handlers("txtprogressbar")  # Text progress bar in the console

# Wrap the call with with_progress to enable progress updates.
# Assume 'coverage_data' is your GRanges object imported from the BigWig.
# Also assume 'repeat_gr' is your GRanges object built from RepeatMasker.
# Here, we process (say) the first 1000 intervals.
metrics_df_vectorized <- compute_interval_metrics_vectorized_progress(
  fasta_file = fasta_path,
  intervals = coverage_data,
  repeat_gr = repeat_gr,
  chunk_size = 10000  # adjust chunk size as needed
)

metrics_df_vectorized <- metrics_df_vectorized[,1:6]

# ggplot(metrics_df_vectorized, aes(x=gc_content))+geom_histogram(bins = 2000)
# ggplot(metrics_df_vectorized, aes(x=repeat_fraction))+geom_histogram(bins = 200)
metrics_df_vectorized$score <- coverage_data$score
# ggplot(metrics_df_vectorized, aes(x=score))+geom_histogram(bins = 10000)
mod <- lm(score~gc_content+repeat_fraction, data = metrics_df_vectorized)
metrics_df_vectorized$predicted_score <- predict(mod, newdata = metrics_df_vectorized)
final <- compute_weighted_log2_ratios(metrics_df_vectorized, observed_col = "score",
                             reference_col = "predicted_score", bin_length_col = "width",
                             repeat_fraction_col = "repeat_fraction")
fds <- final[final$seqnames=="chr9",]
ggplot(fds, aes(x = start, y = log2_ratio))+geom_point(size = 0.1)

metrics_df_vectorized$score<-NULL
gc_repeat_data<-metrics_df_vectorized
usethis::use_data(gc_repeat_data, overwrite = TRUE, compress = "xz") # overwrite = TRUE to overwrite existing files, compress = "xz" for better compression

