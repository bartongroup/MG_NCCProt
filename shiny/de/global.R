# DIFFERENTIAL EXPRESSION EXPLORER
#
# Global settings and variables
#


# set library path for the Shiny server, needs to be stored in the ../libdir.txt
# file - one line of text.

lib_dir <- readLines("../libdir.txt")
if (dir.exists(lib_dir)) .libPaths(lib_dir)

library(shiny)
library(bslib)
library(bsicons)
library(jsonlite)

# datatable needs this option to display Inf correctly

options(htmlwidgets.TOJSON_ARGS = list(na = 'string'))

options(dplyr.summarise.inform = FALSE)
select <- dplyr::select
rename <- dplyr::rename

# Source all R files

files_R <- list.files(c("../R", "../modules"), pattern = "*.R$", full.names = TRUE)
sr_ <- sapply(files_R, source)

# Read config json file

CONFIG <- jsonlite::read_json("../config.json")

# Read all data, make it available to all modules

DATA <- map(names(CONFIG$experiments), function(expmnt) {
  dats <- sh_read_data(c("data", "metadata", "features", "de", "fterms"), expmnt, with_progress = FALSE)
  c(
    dats,
    list(
      id2name = set_names(dats$features$name, dats$features$id),
      name2id = set_names(dats$features$id, dats$features$name)
    )
  )
}) |> 
  set_names(CONFIG$experiments)

EXPERIMENTS <- names(DATA) |> 
  as.character()

