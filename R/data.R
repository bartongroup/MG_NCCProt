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
    mutate(gene_symbol = str_remove(gene_symbols, ";.+$")) |> 
    select(id, gene_symbol)
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


