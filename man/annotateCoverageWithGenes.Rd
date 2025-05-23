% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{annotateCoverageWithGenes}
\alias{annotateCoverageWithGenes}
\title{Annotate Coverage Data with Gene Information}
\usage{
annotateCoverageWithGenes(
  coverage_df,
  species_annotation = org.Hs.eg.db,
  txdb = TxDb.Hsapiens.UCSC.hg38.knownGene,
  return_entrez = FALSE
)
}
\arguments{
\item{coverage_df}{A data frame with at least the following columns: \code{seqnames},
\code{start}, and \code{end} representing coverage intervals.}

\item{txdb}{A transcript database object from which gene information is extracted.
Defaults to \code{TxDb.Hsapiens.UCSC.hg19.knownGene}.}

\item{return_entrez}{Logical. If \code{FALSE} (default), gene symbols are returned; otherwise,
the Entrez gene IDs are returned.}
}
\value{
The original \code{coverage_df} with two additional columns: \code{gene_id} (gene IDs)
  and \code{gene_symbol} (gene symbols or Entrez IDs, depending on \code{return_entrez}).
}
\description{
Converts a coverage data frame into a GRanges object, extracts gene ranges from a transcript
database, finds overlaps between the coverage intervals and genes, and maps gene IDs to gene
symbols (or returns Entrez IDs).
}
\details{
The function converts the input data frame to a \code{GRanges} object using
  \code{GRanges()} and \code{IRanges()}. It then obtains gene ranges via \code{genes()} from the
  transcript database and finds overlaps with \code{findOverlaps()}. Gene IDs are extracted from the
  metadata columns of the gene ranges, and, if requested, mapped to gene symbols using
  \code{mapIds()}.
}
\keyword{internal}
