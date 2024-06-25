library(tidyverse)
library(readxl)
library(limma)
library(dexdash)

source("maxquant.R")
source("de.R")
source("batch.R")


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
) |> 
  remove_batch_effects(formula = "~ group")

# Differential expression

da_full <- limma_de_f(prot, "~ group + batch", what = "abu_med", logfc_limit = 0, fdr_limit = 0.01)
da_pariwise <- limma_de(prot, group_var = "group", what ="abu_limma", logfc_limit = 0, fdr_limit = 0.01)

# Prepare set for dexdash, full modell

dex_full <- dexdash::dexdash_set(
  de = da_full |> 
    select(id, log_fc = logFC, expr = AveExpr, p_value = PValue, contrast = contrast),
  data = prot$dat |> 
    select(id, sample, value = abu_med),
  metadata = des |> 
    select(sample, group, batch, replicate),
  name = "Full model"
)

# Prepare second set for dexdash, pairwise model, data with batch effects removed

dex_pairwise <- dexdash::dexdash_set(
  de = da_pariwise |> 
    select(id, log_fc = logFC, expr = AveExpr, p_value = PValue, contrast = contrast),
  data = prot$dat |> 
    select(id, sample, value = abu_limma),
  metadata = des |> 
    select(sample, group, batch, replicate),
  name = "Pairwise model"
)

# Create a dexdash_list object containing two sets in it:

dex <- dexdash::dexdash_list(dex_full, dex_pairwise)

features <- prot$info |> 
  select(id, name = gene_symbol, description = protein_names)

terms <- download_functional_terms(species = "human")
fterms <- prepare_functional_terms(terms, feature_name = "gene_symbol")
#fterms <- readRDS("fterms.rds")


# Run dexdash, now with new arguments. The first argument is the object
# containing two dexdash sets.

dexdash::run_app(dex, features, fterms)
