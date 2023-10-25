UNIPROT_MAPPING_FILE <- "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/HUMAN_9606_idmapping.dat.gz"
METADATA_FILE <- "info/Design_experiments_4rep.xlsx"

TAXONOMY_ID <- 9606
SPECIES <- "human"
ENSEMBL_DATASET <- "hsapiens_gene_ensembl"
ENSEMBL_VERSION <- "110"

CHROMOSOMES <- c(1:22, "X", "Y", "MT")


CONTRASTS_E1 <- c(
  "TPL_nas-DMSO_nas",
  "DRB_nas-DMSO_nas",
  "TPL_2h-DMSO_2h",
  "DRB_2h-DMSO_2h",
  "TPL_2h-TPL_nas",
  "DRB_2h-DRB_nas",
  "DMSO_2h-DMSO_nas",
  "DMSO_2h-Neg"
)

CONTRASTS_E2 <- c(
  "TPL_6h-DMSO_6h",
  "DRB_6h-DMSO_6h",
  "TPL_2h-DMSO_2h",
  "DRB_2h-DMSO_2h",
  "TPL_6h-TPL_2h",
  "DRB_6h-DRB_2h",
  "DMSO_2h-Neg",
  "DMSO_6h-Neg",
  "TPL_6h-Neg",
  "DRB_6h-Neg"
)

EXPERIMENTS <- tibble::tribble(
  ~name, ~experiment, ~protocol, ~file, ~contrasts,
  "e1_input", "E1", "Input", "mq_data/Experiment 1 Input/proteinGroups_LI.txt", CONTRASTS_E1,
  "e1_ip", "E1", "IP", "mq_data/Exp1 including_4thReplicate/proteinGroups.txt", CONTRASTS_E1,
  "e2_ip", "E2", "IP", "mq_data/Experiment 2/proteinGroups.txt", CONTRASTS_E2
) |> 
  dplyr::mutate(selection = stringr::str_glue("experiment == '{experiment}' & protocol == '{protocol}'"))

TREATMENTS <- tibble::tribble(
  ~name, ~treat, ~ctr, ~ctr_lograt,
  "tpl", "TPL", "treatmentTPL", "TPL_nas-DMSO_nas",
  "drb", "DRB", "treatmentDRB", "DRB_nas-DMSO_nas"
)

PROTEINS_DATA_COLUMNS <- tibble::tribble(
  ~name, ~raw_name, ~type,
  "proteins", "Protein IDs", "c",
  "protein", "Majority protein IDs", "c",
  "gene_symbols", "Gene names", "c",
  "protein_names", "Protein names", "c",
  "sequence_length", "Sequence length", "n",
  "n_razor_unique", "Razor + unique peptides", "n",
  "reverse", "Reverse", "c",
  "contaminant", "Potential contaminant", "c"
)

PROTEINS_ID_COLUMNS <- c("id")
# Using %in% is necessary when the alternative is NA
PROTEINS_FILTER <- "n_razor_unique > 2 & !(reverse %in% '+') & !(contaminant %in% '+')"
MEASURE_COL_PATTERN <- "Reporter intensity corrected \\d{1,2} SB"

FDR_LIMIT <- 0.01
LOGFC_LIMIT <- 0


###-------------------------------------------------------------------

read_metadata <- function(file) {
  d1 <- readxl::read_excel(file, sheet = "Experiment 1") |> 
    rename(sample = `Sample name`, protocol = Sample, replicate = Replicate) |> 
    add_column(experiment = "E1", `6h` = "-") |> 
    mutate(`TMT channel in Maxquant` = as.character(`TMT channel in Maxquant`))
  d2 <- readxl::read_excel(file, sheet = "Experiment 2") |> 
    rename(sample = `Sample name`, protocol = Sample, replicate = Replicate) |>
    mutate(`TMT channel in Maxquant` = as.character(`TMT channel in Maxquant`)) |> 
    add_column(experiment = "E2", nascent = "-")
  d <- bind_rows(d1, d2)
  
  time_cols <- c("Nascent", "1h",  "2h", "6h")
  treatment_cols <- c("DMSO", "TPL", "DRB")
  
  neg <- d |> 
    filter(str_detect(sample, "Neg")) |> 
    select(-all_of(c(time_cols, treatment_cols))) |> 
    mutate(treatment = "Neg", time_point = "Neg")
  meta <- d |> 
    pivot_longer(all_of(time_cols), names_to = "time_point") |>
    filter(value == "+") |> 
    select(-value) |> 
    pivot_longer(all_of(treatment_cols), names_to = "treatment") |>
    filter(value == "+") |> 
    select(-value) 
  bind_rows(neg, meta) |> 
    mutate(time_point = if_else(time_point == "Nascent", "nas", time_point)) |> 
    select(-sample) |> 
    unite(sample, c(experiment, protocol, treatment, time_point, replicate), remove = FALSE) |> 
    mutate(sample = str_replace(sample, "Neg_Neg", "Neg")) |> 
    arrange(experiment, protocol, treatment, time_point) |> 
    clean_names() |> 
    select(experiment, sample, protocol, treatment, time_point, replicate, batch, tmt_channel = tmt_channel_in_maxquant, tmt_tag) |> 
    unite(group, c(treatment, time_point), remove = FALSE) |> 
    mutate(group = str_replace(group, "Neg_Neg", "Neg")) |> 
    mutate(across(c(experiment, sample, protocol, group, treatment, time_point, replicate, batch), as_factor)) |> 
    add_column(bad = FALSE) |> 
    mutate(
      treatment = fct_relevel(treatment, "DMSO"),
      time_point = fct_relevel(time_point, "nas"),
      treatment = fct_relevel(treatment, "Neg", after = Inf)
    ) |> 
    arrange(experiment, protocol, treatment, time_point) |> 
    mutate(group = as_factor(as.character(group)))
}
