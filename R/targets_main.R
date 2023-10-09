targets_main <- function() {
  
  read_metadata <- tar_plan(
    metadata = read_metadata(METADATA_FILE)
  )
  
  annotations <- tar_plan(
    uni_gene = download_uniprot_mapping(UNIPROT_MAPPING_FILE)
  )
  
  
  map_experiments <- tar_map(
    values = EXPERIMENTS,
    names = name,
    
    tar_plan(
      prot = read_mq(file, PROTEINS_DATA_COLUMNS, metadata, uni_gene, sel_meta = selection, filt_data = PROTEINS_FILTER, measure_col_pattern = MEASURE_COL_PATTERN),
      
      fig_detection = plot_detection(prot),
      fig_sample_detection = plot_sample_detection(prot),
      fig_sample_distribution = plot_sample_ridges(prot),
      fig_clustering = plot_clustering(prot, colour_var = "batch"),
      fig_matrix = plot_distance_matrix(prot, min_cor = 0.6),
      fig_pca = plot_pca(prot, shape_var = "batch"),
      
      de_full = limma_de_f(prot, "~ treatment + time_point + batch", what = "abu_med", filt = "treatment != 'Neg'"),
      de_contrasts = limma_de(prot, contrasts = unlist(contrasts)),
      
      fig_volcano_full = plot_volcano(de_full),
      fig_ma_full = plot_ma(de_full),
      fig_pdist_full = plot_pdist(de_full),
      fig_volcano_contrasts = plot_volcano(de_contrasts),
      fig_ma_contrasts = plot_ma(de_contrasts),
      fig_pdist_contrasts = plot_pdist(de_contrasts)
    )
  )
  
  selections <- tar_plan(
    da_genes_contrasts_tpl = de_contrasts_e1_ip |> filter(contrast == "TPL_nas-DMSO_nas" & FDR < 0.05) |> pull(id),
    da_samples_tpl = prot_e1$metadata |> filter(treatment %in% c("DMSO", "TPL")) |> pull(sample),
  )
  
  merged_e1 <- tar_plan(
    prot_e1 = merge_sets(prot_e1_ip, prot_e1_input),
    
    fig_prot_tpl = plot_protein(prot_e1, pids = da_genes_contrasts_tpl, sample_sel = da_samples_tpl)
  )
  
  c(
    annotations,
    read_metadata,
    map_experiments,
    selections,
    merged_e1
  )
}