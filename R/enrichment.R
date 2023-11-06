


fgsea_run <- function(term_data, res, min.size = 3) {
  res <- res |>
    filter(!is.na(value) & !is.na(feature_id)) |> 
    mutate(first_id = str_remove(feature_id, ";.*$"))
  ranks <-  set_names(res$value, res$first_id)
  fgsea::fgsea(pathways = term_data$term2feature, stats = ranks, nproc = 6, minSize = min.size, eps = 0) |>
    as_tibble() |> 
    rename(term_id = pathway) |> 
    mutate(term_name = term_data$term2name[term_id]) |> 
    arrange(NES) |>
    select(term_id, term_name, p_value = pval, fdr = padj, nes = NES, size, leading_edge = leadingEdge)
}

fgsea_groups <- function(d, term_data, feature_var, value_var, group_var) {
  d |>
    mutate(value = get(value_var), feature_id = get(feature_var)) |>
    group_split(!!sym(group_var)) |>
    map_dfr(function(w) {
      fgsea_run(term_data, w) |>
        mutate(!!group_var := dplyr::first(w[[group_var]]))
    })
}

fgsea_all_terms <- function(d, terms, feature_var = "gene_id", value_var = "logFC", group_var = "contrast") {
  ontologies <- names(terms)
  map(ontologies, function(ont) {
    cat(str_glue("  Computing fgsea for {ont}\n\n"))
    fgsea_groups(d, terms[[ont]], feature_var, value_var, group_var)
  }) |>
    set_names(ontologies)
}





plot_fgsea_enrichment <- function(term_id, res, trms, value_var = "logFC") {
  lst <- trms$term2feature[[term_id]]
  rnks <- set_names(res[[value_var]], res$gene_id)
  fgsea::plotEnrichment(lst, rnks)
}

split_genes_fgsea <- function(se, fg, groupvar = "contrast") {
  fg |> 
    filter(padj < 0.05) |> 
    group_split(term, !!sym(groupvar)) |> 
    map_dfr(function(w) {
      term <- as.character(w$term)
      gr <- as.character(w[[groupvar]])
      genes <- w$leading_edge[[1]]
      se |> 
        filter(gene_id %in% genes & !!sym(groupvar) == gr) |> 
        mutate(term_id = term, .before = "gene_id")
    })
}


gsea_de <- function(gse, res, fdr_limit = 0.01) {
  contrasts <- unique(res$contrast) |> 
    as.character()
  ontologies <- names(gse)
  
  map_dfr(ontologies, function(ont) {
    map_dfr(contrasts, function(ctr) {
      res_sel <- res |> 
        filter(contrast == ctr)
      gse_sel <-  gse[[ont]] |> 
        filter(contrast == ctr)
      sig_genes <- res_sel |> 
        filter(sig) |> 
        pull(gene_id)
      
      gse_sel |> 
        filter(fdr < fdr_limit) |> 
        unnest(leading_edge) |> 
        filter(leading_edge %in% sig_genes) |> 
        rename(gene_id = leading_edge) |> 
        left_join(res_sel, by = "gene_id") |> 
        select(term_id, term_name, nes, fdr, gene_id, gene_symbol, logFC, logCPM, PValue, FDR) |> 
        add_column(ontology = ont, .before = 1) |> 
        add_column(contrast = ctr)
    })
  })
}


get_terms_str <- function(gso, query, fdr_limit = 0.05) {
  gso |> 
    filter(str_detect(term_name, query) & padj < fdr_limit) |>
    pull(term_id)
}



gsa_text_selection <- function(gse, texts, fdr_limit = 0.05) {
  ontologies <- names(gse)
  map(ontologies, function(ont) {
    g <- gse[[ont]]
    map(texts, function(txt) {
      g |> 
        filter(str_detect(tolower(term_name), txt) & fdr < fdr_limit) |> 
        add_column(search = txt)
    }) |> 
      list_rbind()
  }) |> 
    list_rbind()
}



gse_example <- function(de, terms, gse) {
  genes_used <- de |>
    pull(gene_id) |>
    unique()
  n_genes_used <- length(genes_used)
  n_genes_in_term <- terms[[GSE_EXAMPLE$ontology]]$mapping |>
    filter(term_id == GSE_EXAMPLE$term_id) |>
    filter(gene_id %in% genes_used) |> 
    nrow()
  n_leading_edge <- gse[[GSE_EXAMPLE$ontology]] |>
    filter(term_id == GSE_EXAMPLE$term_id & contrast == GSE_EXAMPLE$contrast) |>
    pull(leading_edge) |>
    unlist() |>
    length()
  term_name <- terms[[GSE_EXAMPLE$ontology]]$terms |> 
    filter(term_id == GSE_EXAMPLE$term_id) |> 
    pull(term_name)
  list(
    term_id = GSE_EXAMPLE$term_id,
    term_name = term_name,
    contrast = GSE_EXAMPLE$contrast,
    n_genes_used = n_genes_used,
    n_genes_in_term = n_genes_in_term,
    n_leading_edge = n_leading_edge
  )
}