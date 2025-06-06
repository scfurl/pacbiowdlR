% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cna_plot.R
\name{weighted_segment_log2}
\alias{weighted_segment_log2}
\title{Weighted Segment Log2}
\usage{
weighted_segment_log2(segment_data)
}
\arguments{
\item{segment_data}{A \code{data.frame} with columns \code{log2_ratio} and \code{weight}.}
}
\value{
A numeric value representing the weighted average log2 ratio. Returns \code{NA} if the total weight is zero.
}
\description{
Calculates a weighted average log2 ratio for a segment given the per-observation log2 ratios and weights.
}
\details{
The function computes the sum of the weights, and if non-zero, returns the weighted average
  of the log2 ratios. If the total weight is zero, \code{NA} is returned.
}
\examples{
seg_data <- data.frame(log2_ratio = c(0.5, 1.2, -0.7),
                       weight = c(1, 2, 1))
ws <- weighted_segment_log2(seg_data)
print(ws)

}
\keyword{internal}
