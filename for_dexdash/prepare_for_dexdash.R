library(tidyverse)
library(readxl)
library(limma)
library(dexdash)

source("maxquant.R")
source("de.R")


# Read experiment description, format data to match maxquant column names.

des <- read_excel("Description.xlsx") |> 
  select(
    batch = Experiment,
    tmt_channel = Channel,
    group = Samples
  ) |>  
   # Provided Excel file contains `Experiment` names SB40 (channels 1-5) and SB41
   # (channels 6-10), while the MaxQuant data contains only SB40SB41 (channels
   # 1-10). I presume these are the same, so I need to change `batch` variable to
   # match MaxQuant columns
  mutate(batch = if_else(batch %in% c("SB40", "SB41"), "SB40SB41", batch)) |> 
  # Format data to match MaxQuant
  mutate(
    batch = paste0(batch, "-IP"),
    tmt_channel = as.character(tmt_channel),
    group = str_replace(group, "-", "_"),     # minus is a bit awkward in column names
    bad = FALSE  # Used to remove bad data
  ) |> 
  # Need unique sample names
  group_by(group) |> 
  mutate(replicate = row_number()) |> 
  ungroup() |> 
  unite(sample, c(group, replicate), remove = FALSE) |> 
  arrange(group, replicate) |> 
  # Convert `group` to factor, set first level to `Neg_IP` for differential expression
  mutate(group = fct_relevel(group, "Neg_IP"))

# Read protein groups

prot <- read_mq(
  mq_file = "proteinGroups_SB37-41.txt",
  data_cols = PROTEINS_DATA_COLUMNS,
  des = des,
  filt_data = PROTEINS_FILTER,
  measure_col_pattern = MEASURE_COL_PATTERN
)

# Differential expression

da_full <- limma_de_f(prot, "~ group + batch", what = "abu_med", logfc_limit = 0, fdr_limit = 0.01)


# Prepare data for dexdash

de <- da_full |> 
  select(id, log_fc = logFC, expr = AveExpr, p_value = PValue, contrast = contrast)

data <- prot$dat |> 
  select(id, sample, value = abu_med)

metadata <- des |> 
  select(sample, group, batch)

features <- prot$info |> 
  select(id, name = gene_symbol, description = protein_names)

terms <- download_functional_terms(species = "human")
fterms <- prepare_functional_terms(terms, feature_name = "gene_symbol")
#fterms <- readRDS("fterms.rds")


# Run dexdash

run_app(de, data, metadata, features, fterms)
