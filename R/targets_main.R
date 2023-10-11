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
      # read data file
      prot = read_mq(file, PROTEINS_DATA_COLUMNS, metadata, uni_gene, sel_meta = selection, filt_data = PROTEINS_FILTER,
                     measure_col_pattern = MEASURE_COL_PATTERN),
      
      # overview figures
      fig_detection = plot_detection(prot),
      fig_sample_detection = plot_sample_detection(prot),
      fig_sample_distribution = plot_sample_ridges(prot),
      fig_clustering = plot_clustering(prot, colour_var = "batch"),
      fig_matrix = plot_distance_matrix(prot, min_cor = 0.6),
      fig_pca = plot_pca(prot, shape_var = "batch"),
      
      # differential abundance
      da_full = limma_de_f(prot, "~ treatment + time_point + batch", what = "abu_med", filt = "treatment != 'Neg'",
                           logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT),
      da_block = limma_de_block(prot, "~ treatment + time_point", block_var = "batch", filt = "treatment != 'Neg'",
                                logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT),
      da_contrasts = limma_de(prot, contrasts = unlist(contrasts), logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT),
      figs_full = plot_de(da_full),
      figs_block = plot_de(da_block),
      figs_contrasts = plot_de(da_contrasts),
      exp_da_full = export_table(da_full),
      exp_da_contrasts = export_table(da_contrasts),
      
      # per batch log-ratios
      lograt = make_batch_lograt(prot, unlist(contrasts)),
      fig_sample_distribution_lograt = plot_sample_ridges(lograt, what = "logFC", fill_var = "group"),
      fig_clustering_lograt = plot_clustering(lograt, what = "logFC_quant", colour_var = "batch"),
      fig_matrix_lograt = plot_distance_matrix(lograt, min_cor = 0.6, what = "logFC_quant"),
      fig_pca_lograt = plot_pca(lograt,  what = "logFC_quant", colour_var = "group", shape_var = "batch"),
      
      # differential abundance
      dl = limma_de_ratio(lograt, what = "logFC_quant", fdr_limit = FDR_LIMIT),
      figs_dl = plot_de(dl),
      exp_dl = export_table(dl),

      # map through TPL and DRB
      tar_map(
        values = TREATMENTS,
        names = name,
        
        tar_plan(
          da_samples = prot$metadata |> filter(treatment %in% c("DMSO", treat)) |> pull(sample),
          da_pids_full = da_full |> filter(contrast == ctr & sig) |> pull(id),
          da_genes_full = id2gene(da_pids_full, prot$id_prot_gene),
          fig_prots_da_full = plot_protein(prot, pids = da_pids_full, sample_sel = da_samples, ncol = 4),
          
          dl_pids = dl |> filter(contrast == ctr_lograt & sig) |> pull(id),
          fig_prots_dl = plot_protein(prot, pids = dl_pids, sample_sel = da_samples, ncol = 4)
        )
      ),
    )
  )
  
  merged_e1 <- tar_plan(
    prot_e1 = merge_sets(prot_e1_ip, prot_e1_input)
  )
  
  input_normalisation <- tar_plan(
    prot_e1_inpnorm = normalise_to_input(prot_e1_ip, prot_e1_input),
    
    da_full_inpnorm = limma_de_f(prot_e1_inpnorm, "~ treatment + time_point + batch", what = "abu_input", filt = "treatment != 'Neg'",
                                logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT),
    da_contrasts_inpnorm = limma_de(prot_e1_inpnorm, what = "abu_input", contrasts = contrasts_e1_ip, logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT),
    
    figs_full_inpnorm = plot_de(da_full_inpnorm),
    figs_contrasts_inpnorm = plot_de(da_contrasts_inpnorm),
    
    tar_map(
      values = TREATMENTS,
      names = name,
      
      tar_plan(
        da_samples_inpnorm = prot_e1$metadata |> filter(treatment %in% c("DMSO", treat)) |> pull(sample),
        da_pids_full_inpnorm = da_full_inpnorm |> filter(contrast == ctr & sig) |> pull(id),
        fig_prots_da_full_inpnorm = plot_protein(prot_e1_ip, pids = da_pids_full_inpnorm, sample_sel = da_samples_inpnorm)
      )
    )
  )
  
  selections <- tar_plan(
    contrasts_e1_ip = unlist(EXPERIMENTS[2, ]$contrasts),
    da_pids_contrasts_tpl = da_contrasts_e1_ip |> filter(contrast == "TPL_nas-DMSO_nas" & sig) |> pull(id),
    da_pids_contrasts_drb = da_contrasts_e1_ip |> filter(contrast == "DRB_nas-DMSO_nas" & sig) |> pull(id),
    
    da_genes_ip_only_drb = setdiff(da_genes_full_drb_e1_ip, da_genes_full_drb_e1_input),
    da_genes_ip_only_tpl = setdiff(da_genes_full_tpl_e1_ip, da_genes_full_tpl_e1_input)
  )
  
  c(
    annotations,
    read_metadata,
    map_experiments,
    selections,
    merged_e1,
    input_normalisation
  )
}