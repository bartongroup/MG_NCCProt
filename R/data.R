#' Convert a data frame to a matrix
#'
#' This function converts a data frame containing protein abundances to a
#' matrix, with rows representing protein IDs and columns representing samples.
#'
#' @param dat A data frame containing protein abundances
#' @param what A character string specifying which column should be used as
#'   values in the matrix (default: "abu_norm")
#' @param names A character string specifying which column should be used as
#'   column names in the matrix (default: "sample")
#'
#' @return A matrix with protein IDs as row names, sample names as column names,
#'   and the specified values in the cells
dat2mat <- function(dat, what = "abu_norm", names = "sample") {
  dat |> 
    pivot_wider(id_cols = id, names_from = !!names, values_from = !!what) |> 
    column_to_rownames("id") |> 
    as.matrix()
}

mat2dat <- function(tab, to_what = "abu_norm", names = "sample") {
  tab |> 
    as.data.frame() |> 
    rownames_to_column("id") |> 
    pivot_longer(-id, names_to = names, values_to = to_what)
}


add_genes <- function(res, info) {
  g <- info |> 
    select(id, gene_symbols, majority_proteins = protein)
  res |> 
    left_join(g, by = "id")
}

download_uniprot_mapping <- function(uri) {
  cache_file <- "cache/unigene.tsv"
  
  if(!dir.exists("cache"))
    dir.create("cache")
  if(!file.exists(cache_file)) {
    read_tsv(uri, col_names = c("uniprot", "what", "gid"), show_col_types = FALSE) |> 
      filter(what == "Gene_Name") |> 
      select(uniprot, gene_symbol = gid) |> 
      write_tsv(cache_file)
  }
  read_tsv(cache_file, show_col_types = FALSE)
}


id2gene <- function(pids, id_prot_gene) {
  id_prot_gene |>
    filter(id %in% pids) |>
    pull(gene_symbols) |>
    unique()
}



make_batch_lograt <- function(set, contrasts, what = "abu_med") {
  meta <- map(contrasts, function(ctr) {
    groups <- str_split(ctr, "-") |> unlist()
    set$metadata |>
      filter(group %in% groups) |>
      pivot_wider(id_cols = c(batch, experiment, protocol), names_from = group, values_from = sample) |> 
      set_names("batch", "experiment", "protocol", "sample_1", "sample_2") |>
      mutate(across(everything(), as.character)) |> 
      add_column(group = ctr) 
  }) |> 
    list_rbind() |> 
    unite(sample, c(experiment, protocol, group, batch), remove = FALSE) |> 
    mutate(replicate = batch, bad = FALSE)
  
  dat <- map(1:nrow(meta), function(i) {
    r <- meta[i, ]
    set$dat |> 
      filter(sample %in% c(r$sample_1, r$sample_2)) |> 
      mutate(val = get(what)) |> 
      pivot_wider(id_cols = id, names_from = sample, values_from = val) |> 
      drop_na() |> 
      set_names(c("id", "s1", "s2")) |> 
      mutate(logFC = s1 - s2) |> 
      add_column(sample = r$sample, .after = 1) |> 
      select(-c(s1, s2))
  }) |> 
    list_rbind()
  
  # Quantile normalisation
  dat_norm <- dat |>
    pivot_wider(id_cols = id, names_from = sample, values_from = logFC) |> 
    column_to_rownames("id") |> 
    as.matrix() |> 
    preprocessCore::normalize.quantiles(keep.names = TRUE) |> 
    as_tibble(rownames = "id") |> 
    pivot_longer(-id, names_to = "sample", values_to = "logFC_quant") |> 
    mutate(id = as.character(id)) |> 
    drop_na()
  
  dat <- dat |> 
    left_join(dat_norm, by = c("id", "sample"))
  
  list(
    info = set$info,
    id_prot_gene = set$id_prot_gene,
    dat = dat,
    metadata = meta,
    columns = set$columns
  )
}


export_table <- function(df) {
  path <- file.path("tab")
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(df))
  file_name <- file.path(path, str_glue("{obj_name}.csv"))
  df |> 
    mutate(across(where(is.numeric), \(x) {signif(x, 4)})) |> 
    write_csv(file_name)
}


da_table <- function(df) {
  df |> 
    mutate(across(where(is.numeric), \(x) {signif(x, 4)}))
}




remove_batch_effects <- function(set, what = "abu_med", names = "sample",
                      batch_var = "batch", formula = "~ treatment + time_point",
                      filt = "TRUE") {
  
  file_dat <- tempfile("data")
  file_desc <- tempfile("description")

  meta <- set$metadata |> 
    filter(!bad & !!rlang::parse_expr(filt)) |> 
    droplevels() |> 
    mutate(ibatch = as.integer(batch)) |> 
    arrange(ibatch)
  design_mat <- model.matrix(as.formula(formula), data = meta)
  
  dat <- set$dat |> 
    pivot_wider(id_cols = id, names_from = !!names, values_from = !!what)
  dat <- dat[, c("id", as.character(meta$sample))] |> 
    rename(Protein.ID = id)
  tab <- dat |> 
    column_to_rownames("Protein.ID") |> 
    as.matrix()
  
  desc <- meta |> 
    rename(ID = sample) |> 
    mutate(sample = row_number(), batch = ibatch) |> 
    select(ID, sample, batch)
  
  write_tsv(dat, file_dat, na = "NaN")
  write_csv(desc, file_desc)
  
  hr_combat <- HarmonizR::harmonizR(file_dat, file_desc, algorithm = "ComBat")
  hr_limma <- HarmonizR::harmonizR(file_dat, file_desc, algorithm = "limma")
  bat <- limma::removeBatchEffect(tab, batch = meta[[batch_var]], design = design_mat)
  
  dat_bat <- mat2dat(bat, "abu_batch", names)
  dat_combat <- mat2dat(hr_combat, "abu_combat", names)
  dat_limma <- mat2dat(hr_limma, "abu_limma", names)
  set$dat <- set$dat |>   
    inner_join(dat_combat, by = c("id", names)) |> 
    inner_join(dat_limma, by = c("id", names)) |> 
    inner_join(dat_bat, by = c("id", names))
  set$metadata <- meta
  
  set
}


get_ids <- function(da) {
  da |> 
    select(id, gene_symbol = gene_symbols) |> 
    separate_longer_delim(gene_symbol, delim = ";") |> 
    distinct()
}


select_ribosomal_proteins <- function(set, limit = 5) {
  ids <- set$info |> 
    filter(str_detect(protein_names, "ribosomal protein")) |> 
    pull(id)
  
  set$dat |> 
    filter(id %in% ids) |> 
    group_by(id) |> 
    summarise(m = mean(log10(intensity))) |> 
    filter(m > limit) |> 
    pull(id)
}
