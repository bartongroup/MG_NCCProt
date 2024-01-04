biomart_fetch_genes <- function(mart) {
  biomaRt::getBM(attributes = c(
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "gene_biotype",
    "percentage_gene_gc_content",
    "ensembl_gene_id",
    "external_gene_name",
    "entrezgene_id",
    "description"
  ), mart = mart) |>
    dplyr::rename(
      chr = chromosome_name,
      start = start_position,
      end = end_position,
      gene_id = ensembl_gene_id,
      gene_symbol = external_gene_name,
      ncbi_id = entrezgene_id,
      gc_content = percentage_gene_gc_content
    ) |>
    dplyr::mutate(description = str_remove(description, "\\s\\[.*\\]")) |>
    tibble::as_tibble()
}
