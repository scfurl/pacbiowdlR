#' Read and parse JSON file with path substitution for mounted drives
#'
#' @description
#' Reads a JSON file and replaces path prefixes to handle differences between mounted drives
#' across different operating systems or computing environments. This is particularly useful
#' when working with shared files that are accessed through different mount points on
#' different systems (e.g., server paths vs. locally mounted volumes).
#'
#' @param file_path Character; path to the JSON file to be read.
#' @param sub1 Character; the source path prefix to be replaced (default: "/fh/working").
#'   This typically represents the server-side or remote path structure.
#' @param sub2 Character; the target path prefix to replace with (default: "/Volumes").
#'   This typically represents the local mount point or path structure.
#'
#' @return A parsed JSON object with all path strings modified to replace \code{sub1} with \code{sub2}.
#'
#' @details
#' This function handles cross-platform path compatibility issues when JSON files contain
#' absolute paths that need to be translated between different systems. For example, paths that
#' refer to network drives or mounted filesystems might need different prefixes on
#' server environments versus local workstations.
#'
#' The function first checks if the specified file exists, then reads and parses the JSON content.
#' It then recursively processes all string elements in the JSON structure, replacing
#' the specified prefix with an alternative.
#'
#' @examples
#' \dontrun{
#' # Read JSON with server paths and convert to local mounted paths
#' config <- read_json_file("config.json")
#'
#' # Specify custom path substitutions
#' config <- read_json_file("config.json",
#'                         sub1 = "/data/projects",
#'                         sub2 = "/Users/me/mounted/projects")
#' }
#'
#' @importFrom jsonlite fromJSON
#' @export

read_json_file <- function(file_path, sub1 = "/fh/working", sub2 = "/Volumes") {
  # Check if file exists
  if (!file.exists(file_path)) {
    stop("File does not exist: ", file_path)
  }
  # Read and parse JSON
  json_data <- fromJSON(file_path)
  json_data_e <- lapply(json_data, function(loc) gsub(sub1, sub2, loc))
  return(json_data_e)
}
