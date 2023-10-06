read_mq <- function(file, data_cols, id_cols, filt, meta, sel_experiment = NULL) {
  if (!is.null(sel_experiment)) {
    meta <- meta |>
      filter(experiment == sel_experiment)
  }
  
  n2n <- set_names(data_cols$raw_name, data_cols$name)
  raw <- read_tsv(file, col_select = c(data_cols$raw_name, mr$measure_col), show_col_types = FALSE) |>
    rename(all_of(n2n)) |>
    mutate(id = as.character(id))
  dat <- raw |>
    filter(rlang::eval_tidy(rlang::parse_expr(filt))) |>
    select(all_of(id_cols), all_of(mr$measure_col)) |>
    pivot_longer(-all_of(id_cols), names_to = "measure_col", values_to = "value") |>
    mutate(value = na_if(value, 0)) |>
    left_join(mr, by = "measure_col") |>
    select(id, multi, sample, value) |>
    mutate(across(c(id, multi), as.integer)) |>
    drop_na()
  info <- raw |>
    select(-all_of(mr$measure_col)) |>
    mutate(id = as.integer(id))
  
  set <- list(
    info = info,
    dat = dat,
    metadata = meta |> select(-measure_col)
  )
  
  set <- normalise_to_median(set)
  set
}



# At least one non-missing value in all replicates in at least one condition.
get_expressed_ids <- function(set) {
  mr <- set$metadata |>
    select(sample, group) |>
    distinct()
  set$dat |>
    left_join(mr, by = "sample") |>
    group_by(id, group) |>
    summarise(n_tot = n(), n_good = length(na.omit(value))) |>
    ungroup() |>
    group_by(id) |>
    summarise(n_good_conditions = sum(n_tot == n_good)) |>
    filter(n_good_conditions > 0) |>
    select(id)
}


get_info_genes <- function(set) {
  set$info |>
    select(gene_name) |>
    drop_na() |>
    separate_rows(gene_name, sep = ";") |>
    distinct() |>
    pull(gene_name)
}


normalise_to_median <- function(set) {
  med <- set$dat |>
    group_by(sample) |>
    summarise(M = median(value, na.rm = TRUE)) |>
    mutate(M = M / mean(M))
  set$dat <- set$dat |>
    left_join(med, by = "sample") |>
    mutate(value_med = value / M) |>
    select(-M)
  set
}


normalise_to_proteins <- function(pho, pro) {
  # here we ignore a handful of phospho sites that a linked to multiple protein groups
  pho$phospho2prot <- pho$info |>
    select(id, protein_id = protein_ids) |>
    filter(!str_detect(protein_id, ";")) |>
    mutate(protein_id = as.integer(protein_id))
  # add protein intensities and mean protein intensities across conditions to phospho data
  pho$dat <- pho$dat |>
    left_join(pho$phospho2prot, by = "id") |>
    left_join(select(pro$metadata, sample, condition), by = "sample") |>
    left_join(select(pro$dat, protein_id = id, sample, prot = value), by = c("protein_id", "sample")) |>
    group_by(protein_id, condition) |>
    mutate(prot_mean = mean(prot)) |>
    ungroup() |>
    mutate(
      value_prot = value / prot,
      value_prot_mean = value / prot_mean
    ) |>
    select(-c(protein_id, condition, prot, prot_mean))
  pho
}