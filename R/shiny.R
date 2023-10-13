write_shiny_file <- function(subdir, obj, method = c("rds", "qs")) {
  method <- match.arg(method)
  path <- file.path("shiny", "data", subdir)
  if (!dir.exists(path)) dir.create(path)
  obj_name <- deparse(substitute(obj))
  cat(paste("  Writing", obj_name))
  if(method == "qs") {
    file_name <- file.path(path, str_glue("{obj_name}.qs"))
    qs::qsave(obj, file_name)
  } else if(method == "rds") {
    file_name <- file.path(path, str_glue("{obj_name}.rds"))
    write_rds(obj, file_name, compress = "xz")
  }
  cat("\n")
}

save_data_for_shiny <- function(name, dset, de, fterms) {
  gns <- dset$info |> 
    select(id, gene_symbol = gene_symbols, description = protein_names)
  
  de <- de |> 
    select(id, p_value = PValue, fdr = FDR, log_fc = logFC, log_exp = AveExpr, contrast, base) |> 
    left_join(gns, by = "id")
  
  data <- dset$dat
  metadata <- dset$metadata
  features <- dset$info |>
    select(id, name = gene_symbols)
  
  write_shiny_file(name, data)
  write_shiny_file(name, metadata)
  write_shiny_file(name, features)
  write_shiny_file(name, de)
  write_shiny_file(name, fterms)
}




make_star_shiny <- function(star_e1, star_norm) {
  star_e1$dat <- star_e1$dat |> 
    left_join(star_norm$dat |> select(gene_id, sample, logratio_Total))
  star_e1
}
