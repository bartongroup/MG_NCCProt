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



make_batch_lograt <- function(set, contrasts) {
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
      pivot_wider(id_cols = id, names_from = sample, values_from = abu_med) |> 
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
    mutate(id = as.integer(id)) |> 
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

