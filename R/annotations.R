download_uniprot_mapping <- function(uri) {
  read_tsv(uri, col_names = c("uniprot", "what", "gid"), show_col_types = FALSE) |> 
    filter(what == "Gene_Name") |> 
    select(uniprot, gene_symbol = gid)
}


#' Prepare term data for GSEA
#'
#' This is an older version of the function from fenr package, which uses lists,
#' not hashes. We need lists for GSEA.
#'
#' @details
#'
#' Takes two data frames with functional term information (\code{terms}) and
#' gene mapping (\code{mapping}) and converts them into an object required by
#' \code{functional_enrichment} for fast analysis. Terms and mapping can be
#' created with database access functions in this package, for example
#' \code{fetch_reactome} or \code{fetch_go_from_go}.
#'
#' @param terms Information about term names/descriptions. A tibble with columns
#'   \code{term_id} and \code{term_name}.
#' @param mapping Information about term-feature mapping. A tibble with
#'   \code{term_id} and a feature id, as identified with \code{feature_name}
#'   argument. For example, if this tibble contains \code{gene_symbol} and
#'   \code{term_id}, then you need to set \code{feature_name = "gene_symbol"}.
#' @param all_features A vector with all feature ids used as background for
#'   enrichment. If not specified, all features found in \code{mapping} will be
#'   used, resulting in a larger object size.
#' @param feature_name Which column to use from mapping table, e.g.
#'   \code{gene_symbol} or \code{ensembl_gene_id}.
#'
#' @return An object class \code{fenr_terms} required by
#'   \code{functional_enrichment}.
#' @export
#'
#' @examples
#' data(exmpl_all)
#' bp <- fetch_bp()
#' bp_terms <- prepare_for_enrichment(bp$terms, bp$mapping, exmpl_all, feature_name = "gene_symbol")
prepare_for_gsea <- function(terms, mapping, all_features = NULL, feature_name = "gene_id") {
  # Binding variables from non-standard evaluation locally
  feature_id <- term_id <- NULL
  
  # Check terms
  if (!all(c("term_id", "term_name") %in% colnames(terms)))
    stop("Column names in 'terms' should be 'term_id' and 'term_name'.")
  
  # Check mapping
  if (!("term_id" %in% colnames(mapping)))
    stop("'mapping' should contain a column named 'term_id'.")
  
  # Check for feature name
  if (!(feature_name %in% colnames(mapping)))
    stop(feature_name, " column not found in mapping table. Check feature_name argument.")
  
  # Replace empty all_features with everything from mapping
  map_features <- mapping[[feature_name]] |>
    unique()
  if (is.null(all_features)) {
    all_features <- map_features
  } else {
    # Check if mapping is contained in all features
    if (length(intersect(all_features, map_features)) == 0)
      stop("No overlap between 'all_features' and features found in 'mapping'. Did you provide correct 'all_features'?")
  }
  
  # Check for missing term descriptions
  mis_term <- setdiff(mapping$term_id, terms$term_id)
  if (length(mis_term) > 0) {
    dummy <- tibble::tibble(
      term_id = mis_term,
      term_name = rep(NA_character_, length(mis_term))
    )
    terms <- dplyr::bind_rows(terms, dummy)
  }
  
  # List to select term name
  term2name <- terms$term_name |>
    purrr::set_names(terms$term_id)
  
  # feature-term tibble
  feature_term <- mapping |>
    dplyr::rename(feature_id = !!feature_name) |>
    dplyr::filter(feature_id %in% all_features)
  
  # Feature to terms conversion list
  feature2term <- feature_term |>
    dplyr::group_by(feature_id) |>
    dplyr::summarise(terms = list(term_id)) |>
    tibble::deframe()
  
  # Term to feature conversion list
  term2feature <- feature_term |>
    dplyr::group_by(term_id) |>
    dplyr::summarise(features = list(feature_id)) |>
    tibble::deframe()
  
  list(
    term2name = term2name,
    term2feature = term2feature,
    feature2term = feature2term
  ) |>
    structure(class = "fenr_terms")
}


download_functional_terms <- function(species) {
  SPECIES <- list(
    human = list(
      go = "goa_human",
      re = "Homo sapiens",
      kg = "hsa",
      wi = "Homo sapiens"
    ),
    yeast = list(
      go = "sgd",
      re =  "Saccharomyces cerevisiae",
      kg = "sce"
    )
  )
  
  sp <- SPECIES[[species]]
  cat("Loading GO term data\n")
  go <- fenr::fetch_go(species = sp$go)
  cat("Loading Reactome data\n")
  re <- fenr::fetch_reactome(sp$re)
  cat("Loading KEGG data\n")
  kg <- fenr::fetch_kegg(sp$kg)
  #cat("Loading WikiPathways data\n")
  #wi <- fenr::fetch_wiki(sp$wi)
  
  list(
    go = go,
    re = re,
    kg = kg
    #wi = wi
  )
}

# based on protein id, still requires ensembl gene_id -> gene symbol conversion for Reactome
index_functional_terms <- function(terms, sym2id) {
  
  sym2id <- sym2id |> 
    select(id, gene_symbol) |> 
    distinct()
  
  add_protein_id <- function(mapping) {
    mapping |> 
      select(gene_symbol, term_id) |> 
      distinct() |> 
      inner_join(sym2id, by = "gene_symbol", multiple = "all", relationship = "many-to-many") |> 
      select(-gene_symbol) |> 
      distinct()
  }
  
  terms$go$mapping <- terms$go$mapping |> 
    add_protein_id()
  
  terms$re$mapping <- terms$re$mapping |> 
    add_protein_id()
  
  terms$kg$mapping <- terms$kg$mapping |> 
    add_protein_id()
  
  #wi$mapping <- wi$mapping |> 
  #  rename(gene_symbol = text_label) |> 
  #  add_protein_id()
  
  terms
}



prepare_terms_gsea <- function(terms, all_features, feature_name = "id") {
  ontologies <- names(terms)
  
  map(ontologies, function(ont) {
    trm <- terms[[ont]]
    prepare_for_gsea(trm$terms, trm$mapping, all_features, feature_name)
  }) |> 
    set_names(ontologies)
}

prepare_terms_fenr <- function(terms, all_features = NULL, feature_name = "id") {
  ontologies <- names(terms)
  
  map(ontologies, function(ont) {
    trm <- terms[[ont]]
    fenr::prepare_for_enrichment(trm$terms, trm$mapping, all_features, feature_name = feature_name)
  }) |> 
    set_names(ontologies)
}


read_gmt <- function(file_name) {
  assertthat::assert_that(file.exists(file_name))
  
  raw <- read_lines(file_name) |> 
    str_split("\\t")
  
  terms <- raw |> 
    map(function(v) {
      c(
        term_id = v[1],
        term_name = v[2]
      )
    })
  # Old fashioned do.call is much faster than creating individual tibbles in map
  terms <- do.call(rbind, terms) |> 
    as_tibble()
  
  
  mapping <- raw |> 
    map(function(v) {
      c(
        term_id = v[1],
        feature_id = list(v[seq(3, length(v))])
      )
    }) 
  mapping <- do.call(rbind, mapping) |>
    as_tibble() |>
    unnest(c(term_id, feature_id)) |> 
    mutate(across(everything(), as.factor))
  
  list(
    terms = terms,
    mapping = mapping
  )
}
