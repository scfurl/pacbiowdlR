#' @importFrom jsonlite fromJSON
#' @export
#' @keywords internal
#'
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

# Example usage:
# json_content <- read_json_file("path/to/your/file.json")
# print(json_content)
