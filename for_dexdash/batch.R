remove_batch_effects <- function(set, what = "abu_med", names = "sample",
                                 batch_var = "batch", formula = "~ group",
                                 filt = "TRUE") {
  
  file_dat <- tempfile("data")
  file_desc <- tempfile("description")
  
  meta <- set$metadata |> 
    filter(!bad & !!rlang::parse_expr(filt)) |> 
    droplevels() |> 
    mutate(
      batch = as_factor(batch),
      ibatch = as.integer(batch)
    ) |> 
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



mat2dat <- function(tab, to_what = "abu_norm", names = "sample") {
  tab |> 
    as.data.frame() |> 
    rownames_to_column("id") |> 
    pivot_longer(-id, names_to = names, values_to = to_what)
}
