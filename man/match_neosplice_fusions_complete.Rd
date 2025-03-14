% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fusionplot.R
\name{match_neosplice_fusions_complete}
\alias{match_neosplice_fusions_complete}
\title{Match sequenced fusions with NeoSplice database with complete match detection}
\usage{
match_neosplice_fusions_complete(
  seq_input,
  db_file,
  output_file = NULL,
  breakpoint_tolerance = 1e+05,
  sequence_similarity_threshold = 60,
  verbose = TRUE,
  seq_column = "predicted_fusion_sequence",
  max_candidates = 50
)
}
\arguments{
\item{seq_input}{Either a file path to sequenced fusions or a data frame containing fusion data}

\item{db_file}{Path to the NeoSplice database file (can be gzipped)}

\item{output_file}{Optional path to save results as CSV (NULL to skip saving)}

\item{breakpoint_tolerance}{Number of base pairs of tolerance when matching
breakpoint positions (default: 100000)}

\item{sequence_similarity_threshold}{Minimum percent identity required for sequence match (default: 60)}

\item{verbose}{Whether to print progress messages (default TRUE)}

\item{seq_column}{Name of the column containing predicted fusion sequences (default: "predicted_fusion_sequence")}

\item{max_candidates}{Maximum number of candidates to test for sequence similarity (default: 50)}
}
\value{
Data frame of fusion matches with detailed information including complete matches
}
\description{
Enhanced version of match_neosplice_fusions that identifies matches meeting all criteria:
gene names, positions, and sequence similarity.
}
