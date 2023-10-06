targets_main <- function() {
  
  read_metadata <- tar_plan(
    metadata = read_metadata(METADATA_FILE)
  )
  
  read_data <- tar_map(
    values = EXPERIMENTS,
    names = name,
    
    tar_plan(
      prot = read_mq(file, PROTEINS_DATA_COLUMNS, metadata, sel_meta = selection, filt_data = PROTEINS_FILTER, measure_col_pattern = MEASURE_COL_PATTERN),
      
      fig_detection = plot_detection(prot),
      fig_sample_detection = plot_sample_detection(prot),
      fig_sample_distribution = plot_sample_ridges(prot),
      fig_clustering = plot_clustering(prot),
      fig_matrix = plot_distance_matrix(prot),
      fig_pca = plot_pca(prot),
      
      de = limma_de_f(prot_e1_ip, "~ treatment + time_point + batch", what = "abu_med", filt = "treatment != 'Neg'"),
      
      fig_volcano = plot_volcano(de),
      fig_ma = plot_ma(de),
      fig_pdist = plot_pdist(de)
    )
  )
  
  c(
    read_metadata,
    read_data
  )
}