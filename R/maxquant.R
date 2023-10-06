read_mq <- function(file, data_cols, meta, sel_meta, filt_data, measure_col_pattern) {
  meta <- meta |>
    filter(rlang::eval_tidy(rlang::parse_expr(sel_meta)))

  # All measure (reporter) columns
  all_cols <- read_tsv(file, n_max = 0, show_col_types = FALSE) |>
    names() 
  data_cols <- data_cols |> 
    filter(raw_name %in% all_cols)
  
  # We have several tags per TMT reporter, these are replicates
  measure_cols <- tibble(raw_name = all_cols |> str_subset(measure_col_pattern)) |> 
    mutate(tmt_channel = str_extract(raw_name, "(?<=ted )\\d{1,2}")) |> 
    group_by(tmt_channel) |> 
    mutate(replicate = as.character(seq(1, n()))) |> 
    ungroup() |> 
    inner_join(meta, by = c("tmt_channel", "replicate")) |> 
    mutate(type = "n") |> 
    select(name = sample, raw_name, type)
  
  cols <- bind_rows(data_cols, measure_cols)
  # types <- paste0(cols$type, collapse = "")
  
  n2n <- set_names(cols$raw_name, cols$name)
  raw <- read_tsv(file, col_select = cols$raw_name, guess_max = 10000, show_col_types = FALSE) |>
    rename(all_of(n2n))
  # Some files have these missing
  if(all(c("n_razor_unique", "reverse", "contaminant") %in% all_cols)) {
    raw <- raw |> 
      filter(rlang::eval_tidy(rlang::parse_expr(filt_data)))
  }
  raw <- raw |> 
    mutate(id = row_number(), .before = 1)
  
  dat <- raw |>
    select(id, all_of(measure_cols$name)) |>
    pivot_longer(-id, names_to = "sample", values_to = "intensity") |>
    mutate(value = na_if(intensity, 0)) |>
    drop_na()
  info <- raw |>
    select(-all_of(measure_cols$name))

  set <- list(
    info = info,
    dat = dat,
    metadata = meta
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
    summarise(M = median(log10(intensity), na.rm = TRUE)) |>
    mutate(M = M / mean(M))
  set$dat <- set$dat |>
    left_join(med, by = "sample") |>
    mutate(abu_med = log10(intensity) / M) |>
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