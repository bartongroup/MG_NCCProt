targets_main <- function() {
  
  read_data <- tar_plan(
    metadata = read_metadata(METADATA_FILE)
  )
  
  c(
    read_data
  )
}