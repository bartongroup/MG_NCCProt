targets_main <- function() {
  
  read_metadata <- tar_plan(
    metadata = read_metadata(METADATA_FILE)
  )
  
  annotations <- tar_plan(
    mart = useEnsembl(biomart = "ensembl", dataset = ENSEMBL_DATASET, version = ENSEMBL_VERSION),
    terms = download_functional_terms(SPECIES),
    bm_genes = biomart_fetch_genes(mart) |> filter(chr %in% CHROMOSOMES),
    uni_gene = download_uniprot_mapping(UNIPROT_MAPPING_FILE)
  )
  
  
  map_experiments <- tar_map(
    values = EXPERIMENTS,
    names = name,
    
    tar_plan(
      # read data file
      prot_all = read_mq(file, PROTEINS_DATA_COLUMNS, metadata, uni_gene, sel_meta = selection, filt_data = PROTEINS_FILTER,
                     measure_col_pattern = MEASURE_COL_PATTERN),
      prot = remove_batch_effects(prot_all, formula = "~ treatment + time_point"),
      
      # overview figures
      fig_detection = plot_detection(prot),
      fig_sample_detection = plot_sample_detection(prot),
      fig_sample_distribution = plot_sample_ridges(prot, what = "abu_batch"),
      fig_clustering = plot_clustering(prot, colour_var = "treatment", what = "abu_batch"),
      fig_matrix = plot_distance_matrix(prot, min_cor = 0.6, what = "abu_batch"),
      fig_pca = plot_pca(prot, shape_var = "batch", what = "abu_batch") + geom_text_repel(aes(label = time_point)),
      
      fig_clustering_med = plot_clustering(prot_all, colour_var = "batch", what = "abu_med"),
      fig_matrix_med = plot_distance_matrix(prot_all, min_cor = 0.6, what = "abu_med"),
      fig_pca_med = plot_pca(prot, shape_var = "batch", what = "abu_med")  + geom_text_repel(aes(label = time_point)),
      fig_pca_batch = plot_pca(prot, shape_var = "batch", what = "abu_batch")  + geom_text_repel(aes(label = time_point)),
      fig_pca_limma = plot_pca(prot, shape_var = "batch", what = "abu_limma")  + geom_text_repel(aes(label = time_point)),
      fig_pca_combat = plot_pca(prot, shape_var = "batch", what = "abu_combat")  + geom_text_repel(aes(label = time_point)),
      
      
      # differential abundance
      da_full = limma_de_f(prot, "~ treatment + time_point", what = "abu_limma", filt = "treatment != 'Neg'",
                           logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT, base = "Full model"),
      da_int = limma_de_f(prot, "~ treatment * time_point", what = "abu_limma", filt = "treatment != 'Neg'",
                           logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT, base = "Full model"),
      da_contrasts = limma_de(prot, what = "abu_limma", contrasts = unlist(contrasts),
                              logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT, base = "Contrast selection"),

      figs_full = plot_de(da_full),
      figs_contrasts = plot_de(da_contrasts),
      figs_int = plot_de(da_int),
      exp_da_full = export_table(da_full),
      exp_da_contrasts = export_table(da_contrasts),
      
      # interacting proteins
      sel_int = da_int |> filter(sig & str_detect(contrast, ":")),
      fig_prots_int = plot_protein(prot, what = "abu_batch", pids = sel_int, ncol = 4),
      
      # per batch log-ratios
      lograt = make_batch_lograt(prot, what = "abu_batch", unlist(contrasts)),
      fig_sample_distribution_lograt = plot_sample_ridges(lograt, what = "logFC", fill_var = "group"),
      fig_clustering_lograt = plot_clustering(lograt, what = "logFC_quant", colour_var = "batch"),
      fig_matrix_lograt = plot_distance_matrix(lograt, min_cor = 0.6, what = "logFC_quant"),
      fig_pca_lograt = plot_pca(lograt,  what = "logFC_quant", colour_var = "group", shape_var = "batch"),
      
      # differential abundance
      dl = limma_de_ratio(lograt, what = "logFC_quant", fdr_limit = FDR_LIMIT, logfc_limit = LOGFC_LIMIT,
                          base = "Log ratio"),
      figs_dl = plot_de(dl),
      exp_dl = export_table(dl),
      
      # Shiny
      all_ids = prot$info$id,
      iterms = index_functional_terms(terms, bm_genes, prot$id_prot_gene),
      fterms = prepare_terms_fenr(iterms, all_ids),
      da_shiny = bind_rows(da_full, da_contrasts, dl),
      sav_shiny = save_data_for_shiny(name, prot, da_shiny, fterms),
      
      # map through TPL and DRB
      tar_map(
        values = TREATMENTS,
        names = name,
        
        tar_plan(
          da_samples = prot$metadata |> filter(treatment %in% c("Neg", "DMSO", treat)) |> pull(sample),
          da_pids_full = da_full |> filter(contrast == ctr & sig) |> pull(id),
          da_genes_full = id2gene(da_pids_full, prot$id_prot_gene),
          fig_prots_da_full = plot_protein(prot, what = "abu_limma", pids = da_pids_full, sample_sel = da_samples, ncol = 4),
          
          dl_pids = dl |> filter(contrast == ctr_lograt & sig) |> pull(id),
          fig_prots_dl = plot_protein(prot, what = "abu_limma", pids = dl_pids, sample_sel = da_samples, ncol = 4)
        )
      ),
    )
  )
  
  input_normalisation <- tar_plan(
    prot_e1_inpnorm = normalise_to_input(prot_e1_ip, prot_e1_input, what = "abu_batch"),
    
    da_full_inpnorm = limma_de_f(prot_e1_inpnorm, "~ treatment + time_point", what = "abu_input",
                                 logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT, filt = "treatment != 'Neg'",
                                 base = "Full model"),
    da_contrasts_inpnorm = limma_de(prot_e1_inpnorm, what = "abu_input", contrasts = contrasts_e1_ip,
                                    logfc_limit = LOGFC_LIMIT, fdr_limit = FDR_LIMIT, base = "Selected contrasts"),
    
    figs_full_inpnorm = plot_de(da_full_inpnorm),
    figs_contrasts_inpnorm = plot_de(da_contrasts_inpnorm),
    
    # Shiny
    all_ids_inpnorm = prot_e1_inpnorm$info$id,
    iterms_inpnorm = index_functional_terms(terms, bm_genes, prot_e1_inpnorm$id_prot_gene),
    fterms_inpnorm = prepare_terms_fenr(iterms_inpnorm, all_ids_inpnorm),
    da_shiny_inpnorm = bind_rows(da_full_inpnorm, da_contrasts_inpnorm),
    sav_shiny_inpnorm = save_data_for_shiny("e1_inpnorm", prot_e1_inpnorm, da_shiny_inpnorm, fterms_inpnorm),
    
    
    tar_map(
      values = TREATMENTS,
      names = name,
      
      tar_plan(
        da_samples_inpnorm = prot_e1_inpnorm$metadata |> filter(treatment %in% c("DMSO", treat)) |> pull(sample),
        da_pids_full_inpnorm = da_full_inpnorm |> filter(contrast == ctr & sig) |> pull(id),
        fig_prots_da_full_inpnorm = plot_protein(prot_e1_inpnorm, what = "abu_input", pids = da_pids_full_inpnorm,
                                                 sample_sel = da_samples_inpnorm, ncol = 4)
      )
    )
  )
  
  selections <- tar_plan(
    all_ids = prot_e1_input$info$id,
    
    contrasts_e1_ip = unlist(EXPERIMENTS[2, ]$contrasts),
    da_pids_contrasts_tpl = da_contrasts_e1_ip |> filter(contrast == "TPL_nas-DMSO_nas" & sig) |> pull(id),
    da_pids_contrasts_drb = da_contrasts_e1_ip |> filter(contrast == "DRB_nas-DMSO_nas" & sig) |> pull(id),
    da_pids_contrasts_neg = da_contrasts_e1_ip |> filter(contrast == "DMSO_2h-Neg" & sig) |> pull(id),
    
    da_genes_ip_only_drb = setdiff(da_genes_full_drb_e1_ip, da_genes_full_drb_e1_input),
    da_genes_ip_only_tpl = setdiff(da_genes_full_tpl_e1_ip, da_genes_full_tpl_e1_input),
    
    da_ip = da_full_e1_ip |>
      filter(sig & str_detect(contrast, "treatment")) |>
      mutate(
        Input = gene_symbols %in% da_genes_full_drb_e1_input,
        Negative = id %in% da_pids_contrasts_neg
      ) |> 
      left_join(prot_e1_ip$info),
    
    n_ip_all = prot_e1_ip$info |> nrow(),
    n_inpnorm_unique = prot_e1_inpnorm$mapping |> filter(n == 1 & !is.na(id.y)) |> pull(id.x) |> unique() |> length(),
    n_inpnorm_none = prot_e1_inpnorm$mapping |> filter(n == 0) |> pull(id.x) |> unique() |> length(),
    n_inpnorm_multi = prot_e1_inpnorm$mapping |> filter(n > 1) |> pull(id.x) |> unique() |> length()
  )
  
  c(
    annotations,
    read_metadata,
    map_experiments,
    selections,
    input_normalisation
  )
}