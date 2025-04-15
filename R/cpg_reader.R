#' Subset CpG Sites by GRanges Files
#'
#' This function processes a vector of filenames containing CpG site data and subsets them
#' based on overlaps with a provided GRanges object. It returns either a list of raw data.tables
#' from each file or a metric matrix of (averaged) beta values.
#'
#' @param filenames A character vector of filenames to be processed.
#' @param gr A GRanges object. The names of \code{gr} (accessed via \code{names(gr)}) will be used to label
#'   the rows in the metric matrix.
#' @param output A character string indicating the type of output to return. \code{"raw"} (the default)
#'   returns a list of raw data.tables with each table corresponding to a file. \code{"metric"}
#'   returns a matrix (or data.frame) where rows correspond to \code{names(gr)} and columns correspond
#'   to files, with cells containing the (averaged) beta values.
#'
#' @return If \code{output = "raw"}, a list of data.table objects is returned, where each data.table
#'   has an additional column \code{file} indicating the source filename. If \code{output = "metric"},
#'   a matrix (or data.frame) is returned with rownames corresponding to \code{names(gr)} and columns
#'   corresponding to the files. For multiple overlaps in a given GRanges name, the beta values are averaged.
#'
#' @details For each file in \code{filenames}, the function reads in the data using \code{fread()} from the
#'   \code{data.table} package and then subsets the CpG sites that overlap with the GRanges object by calling
#'   \code{subset_cpg_by_GR}. Files yielding no overlaps are removed from the final results. For \code{"metric"}
#'   output, the function creates a vector (with \code{NA} for missing overlaps) for each file, computing the
#'   mean beta value for cases of multiple overlaps.
#'
#' @importFrom data.table fread
#' @importFrom GenomicRanges GRanges
#' @keywords internal
#' @examples
#' \dontrun{
#' # Create an example GRanges object.
#' gr <- GRanges(seqnames = "chr1", ranges = IRanges(start = c(1, 101), width = 100),
#'               names = c("probe1", "probe2"))
#'
#' # Example vector of file names.
#' file_vec <- c("file1.csv", "file2.csv")
#'
#' # Get raw output (list of data.tables)
#' raw_output <- subset_cpg_by_GR_files(file_vec, gr, output = "raw")
#'
#' # Get metric matrix output
#' metric_output <- subset_cpg_by_GR_files(file_vec, gr, output = "metric")
#' }
#'
#' @export
subset_cpg_by_GR_files <- function(filenames, gr, output = c("raw", "metric")) {
  output <- match.arg(output)

  # Process each file and store the raw output in a list.
  raw_list <- lapply(filenames, function(fn) {
    message("Processing file: ", fn)
    dt <- fread(fn)
    res <- subset_cpg_by_GR(dt, gr)
    if (!is.null(res)) {
      # Add a column indicating the source filename.
      res[, file := fn]
    }
    return(res)
  })

  # Remove any NULL results (in case some files had no overlaps).
  names(raw_list) <- filenames
  raw_list <- Filter(Negate(is.null), raw_list)

  if (output == "raw") {
    return(raw_list)
  } else if (output == "metric") {
    # Create a metric matrix.
    # Rows: names from the GRanges object.
    # Columns: filenames; each cell will be the (averaged) beta value.
    probes <- names(gr)
    # Create an empty list to hold beta values for each file.
    metric_list <- vector("list", length(filenames))
    names(metric_list) <- filenames

    for (fn in filenames) {
      dt <- raw_list[[fn]]
      # Initialize a vector of NA with names equal to probes.
      beta_vec <- setNames(rep(NA, length(probes)), probes)
      if (!is.null(dt)) {
        # In case multiple overlaps occur for the same GR name, compute the average beta.
        dt_avg <- dt[, .(beta_val = mean(beta, na.rm = TRUE)), by = name]
        # Update the vector with available beta values.
        beta_vec[dt_avg$name] <- dt_avg$beta_val
      }
      metric_list[[fn]] <- beta_vec
    }
    # Combine the vectors into a matrix.
    metric_mat <- do.call(cbind, metric_list)
    return(metric_mat)
  }
}

#' Recursively Compute the Intersection of Character Vectors
#'
#' This function computes the intersection of a list of character vectors recursively.
#'
#' @param charList A list of character vectors.
#'
#' @return A character vector containing the common elements present in all vectors in \code{charList}.
#'   If \code{charList} is empty, an empty character vector is returned.
#'
#' @details The function handles the following cases:
#'   \itemize{
#'     \item If the list is empty, it returns an empty character vector.
#'     \item If the list contains only one vector, that vector is returned.
#'     \item For multiple vectors, it first computes the intersection of the first two vectors,
#'           then recursively computes the intersection of the result with the rest of the list.
#'   }
#'
#' @examples
#' # Example usage:
#' vectors_list <- list(c("A", "B", "C"), c("B", "C", "D"), c("C", "D", "E"))
#' common_elements <- recursiveIntersect(vectors_list)
#' print(common_elements)  # Expected output might be "C" if that is the only common element.
#' @keywords internal
#' @export
recursiveIntersect <- function(charList) {
  # Base case: If the list is empty, return an empty character vector
  if (length(charList) == 0) {
    return(character(0))
  }

  # Base case: If there is only one vector in the list, return it as the intersection
  if (length(charList) == 1) {
    return(charList[[1]])
  }

  # Recursively intersect the first two vectors.
  currentIntersect <- intersect(charList[[1]], charList[[2]])

  # If only two vectors were provided, currentIntersect is our answer.
  # Otherwise, combine currentIntersect as the first element with the remainder
  # of the list and recursively call the function.
  if (length(charList) == 2) {
    return(currentIntersect)
  } else {
    newList <- c(list(currentIntersect), charList[-c(1,2)])
    return(recursiveIntersect(newList))
  }
}
