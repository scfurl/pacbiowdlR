% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fusionplot.R
\name{get_fusions_from_vcf}
\alias{get_fusions_from_vcf}
\title{Extract Fusions from VCF and Generate a Circos Plot}
\usage{
get_fusions_from_vcf(
  vcf_file,
  gtf_file = NULL,
  genome = "hg38",
  title = "",
  thresh = 20,
  allowed = paste0("chr", c(1:22)),
  highlight = NULL,
  plot = FALSE,
  filter = c("pass", "fully_spanned", "has_strand"),
  annotate = FALSE,
  verbose = F,
  tss_upstream = 2000,
  tss_downstream = 200
)
}
\arguments{
\item{vcf_file}{Character. File path to the gzipped VCF file.}

\item{gtf_file}{Character or \code{NULL}. Optional file path to a GTF file for breakpoint annotation. Required when \code{annotate} is \code{TRUE}.}

\item{genome}{Character. Genome assembly version (default: \code{"hg38"}).}

\item{title}{Character. Title for the plot when \code{plot = TRUE} (default: an empty string).}

\item{thresh}{Numeric. Threshold for read depth filtering (default: \code{20}).}

\item{allowed}{Character vector. Allowed chromosomes (default: \code{paste0("chr", 1:22)}).}

\item{highlight}{Optional. Indices or logical vector of fusion events to highlight in the plot (default: \code{NULL}).}

\item{plot}{Logical. If \code{TRUE} a Circos plot will be generated (default: \code{FALSE}).}

\item{filter}{Character vector. Filtering criteria; supported values include \code{"pass"}, \code{"fully_spanned"}, \code{"protein_coding"}, and \code{"has_strand"}.}

\item{annotate}{Logical. If \code{TRUE}, fusion breakpoints will be annotated using gene coordinates (default: \code{FALSE}).}

\item{verbose}{Logical. if \code{TRUE}, will provide messaging.}

\item{tss_upstream}{Numeric. Number of base pairs upstream of the transcription start site (default: \code{2000}).}

\item{tss_downstream}{Numeric. Number of base pairs downstream of the transcription start site (default: \code{200}).}
}
\value{
A data frame containing fusion events with their breakpoint coordinates and (if requested) annotation details.
}
\description{
Reads a gzipped VCF file containing structural variant data, applies several filters,
optionally annotates fusion breakpoints using a GTF file, and (optionally) generates a Circos plot.
}
