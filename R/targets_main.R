targets_main <- function() {
  
  read_metadata <- tar_plan(
    metadata = read_metadata(METADATA_FILE)
  )
  
  read_data <- tar_map(
    values = EXPERIMENTS,
    names = name,
    
    tar_plan(
      prot = read_mq(file, PROTEINS_DATA_COLUMNS, metadata, sel_meta = selection, filt_data = PROTEINS_FILTER, measure_col_pattern = MEASURE_COL_PATTERN)
    )
  )
  
  c(
    read_metadata,
    read_data
  )
}