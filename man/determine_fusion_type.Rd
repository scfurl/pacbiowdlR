% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fusionplot2.R
\name{determine_fusion_type}
\alias{determine_fusion_type}
\title{Determine the type of gene fusion}
\usage{
determine_fusion_type(
  contig1,
  contig2,
  direction1,
  direction2,
  breakpoint1,
  breakpoint2
)
}
\arguments{
\item{contig1}{Chromosome of gene1}

\item{contig2}{Chromosome of gene2}

\item{direction1}{Direction of gene1 in fusion}

\item{direction2}{Direction of gene2 in fusion}

\item{breakpoint1}{Genomic coordinate of breakpoint in gene1}

\item{breakpoint2}{Genomic coordinate of breakpoint in gene2}
}
\value{
Character string describing the fusion type
}
\description{
Determine the type of gene fusion
}
