read_mq <- function(file, data_cols, metadata, uni_gene, sel_meta, filt_data, measure_col_pattern) {
  meta <- metadata |>
    filter(rlang::eval_tidy(rlang::parse_expr(sel_meta)))

  # All measure (reporter) columns
  all_cols <- read_tsv(file, n_max = 0, show_col_types = FALSE) |>
    names() 
  dat_cols <- data_cols |> 
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
  
  cols <- bind_rows(dat_cols, measure_cols)
  # types <- paste0(cols$type, collapse = "")
  
  n2n <- set_names(cols$raw_name, cols$name)
  raw <- read_tsv(file, col_select = cols$raw_name, guess_max = 10000, show_col_types = FALSE) |>
    rename(all_of(n2n))
  # Some files have these missing
  if(all(c("n_razor_unique", "reverse", "contaminant") %in% cols$name)) {
    raw <- raw |> 
      filter(rlang::eval_tidy(rlang::parse_expr(filt_data)))
  }
  raw <- raw |> 
    mutate(id = row_number(), .before = 1)
  
  dat <- raw |>
    select(id, all_of(measure_cols$name)) |>
    pivot_longer(-id, names_to = "sample", values_to = "intensity") |>
    filter(intensity > 0)
  
  info <- raw |>
    select(-all_of(measure_cols$name)) |> 
    mutate(gene_symbols = if_else(is.na(gene_symbols), proteins, gene_symbols))
  
  # gene_symbols is from maxQuant table
  # gene_symbol is from Uniprot mapping of individual proteins
  id_prot_gene <- info |>
    select(id, uniprot = proteins, gene_symbols) |>
    separate_longer_delim(uniprot, delim = ";") |>
    mutate(uniprot = str_remove(uniprot, "\\-\\d+$")) |>
    left_join(uni_gene, by = c("uniprot"), multiple = "all") |>
    mutate(
      gene_symbol = if_else(is.na(gene_symbol), uniprot, gene_symbol)
    )
  

  set <- list(
    info = info,
    id_prot_gene = id_prot_gene,
    dat = dat,
    metadata = meta,
    columns = cols
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
    select(gene_symbols) |>
    drop_na() |>
    separate_rows(gene_symbols, sep = ";") |>
    distinct() |>
    pull(gene_symbols)
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


normalise_to_input <- function(ip, inp, what = "abu_med") {
  s2g_ip <- ip$metadata |> 
    select(sample, group)
  s2g_inp <- inp$metadata |> 
    select(sample, group)
  
  minp <- inp$dat |> 
    mutate(val = get(what)) |> 
    left_join(s2g_inp, by = "sample") |> 
    group_by(id, group) |> 
    summarise(mean_input = mean(val))
  
  ip$dat <- ip$dat |> 
    mutate(val = get(what)) |>
    left_join(s2g_ip, by = "sample") |> 
    left_join(minp, by = c("id", "group")) |> 
    mutate(abu_input = val - mean_input) |> 
    select(-c(val, group, mean_input))
  
  ip
}


merge_sets <- function(set1, set2) {
  sel_ <- function(info) {
    info |> 
      select(id, x = gene_symbols) |> 
      separate_longer_delim(x, delim = ";")
  }
  
  uni1 <- sel_(set1$info)
  uni2 <- sel_(set2$info)
  
  # id mapping
  mp <- left_join(uni1, uni2, by = "x", relationship = "many-to-many") |> 
    select(-x) |> 
    distinct() |> 
    group_by(id.x) |> 
    mutate(n = n()) |>
    ungroup()
  
  dat2 <- set2$dat |> 
    right_join(mp |> filter(n == 1), by = c("id" = "id.y"), relationship = "many-to-many") |> 
    select(-c(id, n)) |> 
    rename(id = id.x) |> 
    relocate(id, .before = 1) |> 
    drop_na()
  
  meta <- bind_rows(
    set1$metadata, 
    set2$metadata
  )
  cols <- bind_rows(set1$columns, set2$columns)
  
  dat <- bind_rows(set1$dat, dat2) |> 
    select(-starts_with("abu"))
  
  set <- list(
    columns = cols,
    info = set1$info,
    id_prot_gene = set1$id_prot_gene,
    dat = dat,
    metadata = meta,
    mapping = mp
  )
  
  set |> normalise_to_median()
}

