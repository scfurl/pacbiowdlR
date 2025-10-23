#' Subset CpG Sites by GRanges Object
#'
#' @param dt A data.table containing CpG site data
#' @param gr A GRanges object defining regions of interest
#' @param fill_missing Logical. If TRUE, include all probes from gr with NA for uncovered probes
#'
#' @return A data.table with overlapping CpG sites (and NAs if fill_missing=TRUE)
#'
#' @export
subset_cpg_by_GR <- function(dt, gr, fill_missing = FALSE) {
  # Convert data.table to GRanges
  dt_gr <- GRanges(
    seqnames = dt$V1,
    ranges = IRanges(start = dt$V2, end = dt$V3)
  )
  
  # Find overlaps
  overlaps <- findOverlaps(dt_gr, gr)
  
  # If no overlaps and not filling missing, return NULL
  if (length(overlaps) == 0 && !fill_missing) {
    return(NULL)
  }
  
  if (length(overlaps) > 0) {
    # Subset the data.table and add the GRanges names and beta values
    result <- dt[queryHits(overlaps), ]
    result[, name := names(gr)[subjectHits(overlaps)]]
    result[, beta := V4]
  } else {
    result <- NULL
  }
  
  # If fill_missing, create complete set with NAs for missing probes
  if (fill_missing) {
    all_probes <- data.table(name = names(gr))
    
    if (!is.null(result)) {
      # Merge to include all probes
      result <- merge(all_probes, result, by = "name", all.x = TRUE)
    } else {
      # No overlaps found, create all NA result
      result <- all_probes
      result[, beta := NA_real_]
    }
  }
  
  return(result)
}

#' Subset CpG Sites by GRanges Files
#'
#' @param filenames A character vector of filenames to be processed
#' @param gr A GRanges object
#' @param output Output type: "raw" or "metric"
#' @param fill_missing Logical. If TRUE, include all probes with NA for uncovered probes
#'
#' @export
subset_cpg_by_GR_files <- function(filenames, gr, output = c("raw", "metric"), fill_missing = TRUE) {
  output <- match.arg(output)

  # Process each file and store the raw output in a list
  raw_list <- lapply(filenames, function(fn) {
    message("Processing file: ", fn)
    dt <- fread(fn)
    res <- subset_cpg_by_GR(dt, gr, fill_missing = fill_missing)
    if (!is.null(res)) {
      # Add a column indicating the source filename
      res[, file := fn]
    }
    return(res)
  })

  # Remove any NULL results (in case some files had no overlaps)
  names(raw_list) <- filenames
  raw_list <- Filter(Negate(is.null), raw_list)

  if (output == "raw") {
    return(raw_list)
  } else if (output == "metric") {
    # Create a metric matrix
    probes <- names(gr)
    metric_list <- vector("list", length(filenames))
    names(metric_list) <- filenames

    for (fn in filenames) {
      dt <- raw_list[[fn]]
      # Initialize a vector of NA with names equal to probes
      beta_vec <- setNames(rep(NA, length(probes)), probes)
      if (!is.null(dt)) {
        # In case multiple overlaps occur for the same GR name, compute the average beta
        dt_avg <- dt[, .(beta_val = mean(beta, na.rm = TRUE)), by = name]
        # Update the vector with available beta values
        beta_vec[dt_avg$name] <- dt_avg$beta_val
      }
      metric_list[[fn]] <- beta_vec
    }
    # Combine the vectors into a matrix
    metric_mat <- do.call(cbind, metric_list)
    
    # If fill_missing is FALSE, remove rows that are all NA
    if (!fill_missing) {
      rows_with_data <- apply(metric_mat, 1, function(x) !all(is.na(x)))
      metric_mat <- metric_mat[rows_with_data, , drop = FALSE]
    }
    
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
