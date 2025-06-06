% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cpg_reader.R
\name{subset_cpg_by_GR_files}
\alias{subset_cpg_by_GR_files}
\title{Subset CpG Sites by GRanges Files}
\usage{
subset_cpg_by_GR_files(filenames, gr, output = c("raw", "metric"))
}
\arguments{
\item{filenames}{A character vector of filenames to be processed.}

\item{gr}{A GRanges object. The names of \code{gr} (accessed via \code{names(gr)}) will be used to label
the rows in the metric matrix.}

\item{output}{A character string indicating the type of output to return. \code{"raw"} (the default)
returns a list of raw data.tables with each table corresponding to a file. \code{"metric"}
returns a matrix (or data.frame) where rows correspond to \code{names(gr)} and columns correspond
to files, with cells containing the (averaged) beta values.}
}
\value{
If \code{output = "raw"}, a list of data.table objects is returned, where each data.table
  has an additional column \code{file} indicating the source filename. If \code{output = "metric"},
  a matrix (or data.frame) is returned with rownames corresponding to \code{names(gr)} and columns
  corresponding to the files. For multiple overlaps in a given GRanges name, the beta values are averaged.
}
\description{
This function processes a vector of filenames containing CpG site data and subsets them
based on overlaps with a provided GRanges object. It returns either a list of raw data.tables
from each file or a metric matrix of (averaged) beta values.
}
\details{
For each file in \code{filenames}, the function reads in the data using \code{fread()} from the
  \code{data.table} package and then subsets the CpG sites that overlap with the GRanges object by calling
  \code{subset_cpg_by_GR}. Files yielding no overlaps are removed from the final results. For \code{"metric"}
  output, the function creates a vector (with \code{NA} for missing overlaps) for each file, computing the
  mean beta value for cases of multiple overlaps.
}
\examples{
\dontrun{
# Create an example GRanges object.
gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 101), width = 100),
              names = c("probe1", "probe2"))

# Example vector of file names.
file_vec <- c("file1.csv", "file2.csv")

# Get raw output (list of data.tables)
raw_output <- subset_cpg_by_GR_files(file_vec, gr, output = "raw")

# Get metric matrix output
metric_output <- subset_cpg_by_GR_files(file_vec, gr, output = "metric")
}

}
\keyword{internal}
