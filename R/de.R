# Convert raw output from limma into a useful format
tabulate_de <- function(fit) {
  coefs <- colnames(fit$coefficients) |> 
    str_subset("Intercept", negate = TRUE)
  map_dfr(coefs, function(ctr) {
    limma::topTable(fit, coef = ctr, number = 1e6, sort.by = "none") |>
      as_tibble(rownames = "id") |>
      add_column(contrast = ctr)
  }) |> 
    drop_na() |> 
    select(-c(t, B)) |>
    rename(FDR = adj.P.Val, PValue = P.Value) |>
    mutate(
      id = as.integer(id),
      logFC = logFC / log10(2),    # convert into log2
      contrast = factor(contrast, levels = coefs)
    )
}

# DE for selected contrasts; if not specified, all pairs of contrasts 
limma_de <- function(set, contrasts = NULL, group_var = "treatment", what = "abu_norm",
                     filt = "TRUE", names = "sample") {
  meta <- set$metadata |> 
    filter(!bad & !!rlang::parse_expr(filt)) |> 
    mutate(group = get(group_var)) |> 
    filter(!is.na(group)) |> 
    droplevels()
  groups <- unique(as.character(meta$group)) |> 
    janitor::make_clean_names()
  design_mat <- model.matrix(~ 0 + group, data = meta)
  colnames(design_mat) <- groups
  
  tab <- dat2mat(set$dat, what, names)[, as.character(meta[[names]])]
  
  if (is.null(contrasts)) {
    contrasts <- expand_grid(x = as_factor(groups), y = as_factor(groups)) |>
      filter(as.integer(x) < as.integer(y)) |>
      unite(contrast, c(y, x), sep = "-") |>
      pull(contrast)
  }
  contrast_mat <- limma::makeContrasts(contrasts = contrasts, levels = design_mat)

  fit <- tab |>
    limma::lmFit(design_mat) |>
    limma::contrasts.fit(contrasts = contrast_mat) |>
    limma::eBayes()
  
  tabulate_de(fit) |> 
    add_genes(set$info)
}



# DE with formula
limma_de_f <- function(set, formula, what = "abu_norm", filt = "TRUE", names = "sample") {
  meta <- set$metadata |> 
    filter(!bad & !!rlang::parse_expr(filt)) |> 
    droplevels()
  
  tab <- dat2mat(set$dat, what, names)[, as.character(meta[[names]])]
  design_mat <- model.matrix(as.formula(formula), data = meta)

  fit <- tab |>
    limma::lmFit(design_mat) |>
    limma::eBayes()
  
  tabulate_de(fit) |> 
    add_genes(set$info)
}



# one-sample limma against zero
limma_de_ratio <- function(df, what = "logFC", id_var = "participant_id", filt = "TRUE") {
  meta <- df$metadata |> 
    filter(!bad & !!rlang::parse_expr(filt)) |> 
    droplevels()
  
  tab <- dat2mat(df$dat, what = what, names = id_var)
  tab <- tab[, as.character(meta[[id_var]])]
  
  design_mat <- cbind(Intercept = rep(1, ncol(tab)))
  fit <- tab |>
    limma::lmFit(design_mat) |>
    limma::eBayes()
  
  res <- limma::topTable(fit, number = 1e6, sort.by = "none") |>
      as_tibble(rownames = "id") |>
      mutate(id = as.integer(id)) |> 
      rename(FDR = adj.P.Val, PValue = P.Value) |>
      select(-c(t, B)) |> 
      add_column(contrast = "ratio") |> 
      drop_na() |> 
      add_genes(df$info)
  
  res
}



de_list <- function(res, group_var, fdr = "FDR", logfc = "logFC", logfc_limit = 0, fdr_limit = 0.05, name = NULL, split_up_down = FALSE) {
  d <- res |>
    filter(!!sym(fdr) < fdr_limit & abs(!!sym(logfc)) >= logfc_limit) |>
    mutate(group = !!sym(group_var))
  if (split_up_down) {
    d <- d |>
      mutate(direction = if_else(!!sym(logfc) > 0, "up", "down")) |>
      unite(group, c(group, direction), sep = ":")
  }
  d <- d |>
    select(group, id) |>
    group_by(group)
  kname <- ifelse(is.null(name), "", paste0(name, ":"))
  ks <- paste0(kname, group_keys(d)[[1]])
  d |>
    distinct() |>
    group_map(~pull(.x, id)) |>
    set_names(ks)
}


pull_proteins <- function(des) {
  des |>
    select(protein, gene_symbol) |>
    distinct() |>
    arrange(gene_symbol) |>
    filter(!is.na(gene_symbol))
}

make_de_genes <- function(de, fdr_limit = 0.05, logfc_limit = 1) {
  list(
    up = pull_proteins(de |> filter(FDR < fdr_limit & logFC >= logfc_limit)),
    down = pull_proteins(de |> filter(FDR < fdr_limit & logFC <= -logfc_limit))
  )
}

make_de_table <- function(de, info) {
  de |>
    select(id, logFC, PValue, FDR) |>
    left_join(info, by = "id") |>
    select(id, logFC, FDR, protein_accessions, gene_symbol) |>
    mutate(across(c(logFC, FDR), ~signif(.x, 3)))
}
