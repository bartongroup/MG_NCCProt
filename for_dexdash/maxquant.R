#' Data columns for protein data
#'
#' @format A tibble with columns:
#' \describe{
#'   \item{name}{Name of the column in the processed data}
#'   \item{raw_name}{Name of the column in the raw data file}
#'   \item{type}{Type of data ('c' for character, 'n' for numeric)}
#' }
PROTEINS_DATA_COLUMNS <- tibble::tribble(
  ~name, ~raw_name, ~type,
  "proteins", "Protein IDs", "c",
  "protein", "Majority protein IDs", "c",
  "gene_symbols", "Gene names", "c",
  "protein_names", "Protein names", "c",
  "sequence_length", "Sequence length", "n",
  "n_razor_unique", "Razor + unique peptides", "n",
  "ibaq", "iBAQ", "n",
  "reverse", "Reverse", "c",
  "contaminant", "Potential contaminant", "c",
  "only_site", "Only identified by site", "c"
)

#' Filter criteria for proteins data
#'
#' A character string specifying the filter criteria for the protein data.
PROTEINS_FILTER <- "n_razor_unique > 2 & !(reverse %in% '+') & !(contaminant %in% '+') & !(only_site %in% '+')"

#' Measurement column pattern
#'
#' A regular expression pattern for identifying measurement columns.
MEASURE_COL_PATTERN <- "Reporter intensity corrected \\d{1,2} SB"

#' Read MaxQuant data
#'
#' This function reads protein groups data from a MaxQuant output file and
#' processes it according to the specified parameters.
#'
#' @param mq_file A string specifying the path to the MaxQuant proteinGroups.txt
#'   file.
#' @param data_cols A tibble specifying the data columns to be used.
#' @param des A tibble containing the experimental design information.
#' @param filt_data A string specifying the filter criteria to be applied to the
#'   data.
#' @param measure_col_pattern A regular expression pattern for identifying
#'   measurement columns.
#'
#' @return A list containing the processed data, including: \item{info}{A tibble
#'   with useful information about proteins and genes.} \item{id_prot_gene}{A
#'   tibble with individual protein IDs and gene symbols.} \item{dat}{A
#'   long-format tibble with sample intensities.} \item{description}{The
#'   experimental design tibble.} \item{columns}{A tibble with column names and
#'   types.}
read_mq <- function(mq_file, data_cols, des, filt_data, measure_col_pattern) {

  # Find all data_cols that are present in the protein groups file
  all_cols <- read_tsv(mq_file, n_max = 0, show_col_types = FALSE) |>
    names() 
  dat_cols <- data_cols |> 
    filter(raw_name %in% all_cols)
  
  # We have several tags per TMT reporter, these are replicates
  # Create a tibble to match measure column names with samples
  measure_cols <- tibble(raw_name = all_cols |> str_subset(measure_col_pattern)) |> 
    mutate(
      batch = str_extract(raw_name, "SB.+$"),
      tmt_channel = str_extract(raw_name, "(?<=ted )\\d{1,2}")
    ) |> 
    inner_join(des, by = c("tmt_channel", "batch")) |> 
    mutate(type = "n") |> 
    select(name = sample, raw_name, type)
  
  # All column we need to read
  cols <- bind_rows(dat_cols, measure_cols)
  
  # Translate raw column name to a new name used here
  n2n <- set_names(cols$raw_name, cols$name)
  
  # Read protein groups data. Read only columns specified by `cols`. Rename them
  # to useful names.
  raw <- read_tsv(mq_file, col_select = cols$raw_name, guess_max = 10000, show_col_types = FALSE) |>
    rename(all_of(n2n))
  # Some files have these missing
  if(all(c("n_razor_unique", "reverse", "contaminant", "only_site") %in% cols$name)) {
    raw <- raw |> 
      filter(rlang::eval_tidy(rlang::parse_expr(filt_data)))  # filter data
  } else {
    raw <- raw |> 
      filter(!str_detect(protein, "^(CON__|REV__)"))
  }
  # Create a unique id based on gene symbol
  raw <- raw |> 
    mutate(
      gene_symbol = if_else(is.na(gene_symbols), protein, gene_symbols),
      gene_symbol = str_remove(gene_symbol, ";.+$"),
    ) |> 
    mutate(id = gene_symbol |> make.unique(), .before = 1)

  # Pivot data to a long format
  dat <- raw |>
    select(id, all_of(measure_cols$name)) |>
    pivot_longer(-id, names_to = "sample", values_to = "intensity") |>
    filter(intensity > 0)
  
  # Useful information about proteins and genes
  info <- raw |>
    select(-all_of(measure_cols$name), gene_symbol)
  
  set <- list(
    info = info,
    dat = dat,
    metadata = des,
    columns = cols
  )
  
  set <- normalise_to_median(set)
  set
}


#' Normalize data to median
#'
#' This function normalizes the intensity data in a dataset to the median intensity.
#'
#' @param set A list containing the data to be normalized. The list should include:
#' \item{dat}{A tibble with sample intensities.}
#'
#' @return The input list with normalized intensity data.
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

