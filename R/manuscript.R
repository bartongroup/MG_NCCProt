add_axis_groups <- function(d, label_var, group_var, gap = 0.1) {
  d <- d |> 
    ungroup() |> 
    mutate(
      label = get(label_var),
      group = get(group_var)
    )    
  
  ag <- d |> 
    select(label, group) |> 
    distinct() |> 
    group_by(group) |> 
    mutate(igroup = cur_group_id()) |> 
    ungroup() |> 
    mutate(x = row_number() + gap * igroup)
  d |> 
    left_join(ag, by = c("group", "label"))
}


mn_plot_feature_heatmap <- function(set, genes, what = "abu_limma", group_mean = FALSE, text_size = 10,
                                    max_z = NULL) {
  
  prots <- set$id_prot_gene |> 
    filter(gene_symbol %in% genes) |> 
    select(id, gene_symbol) |> 
    distinct()
  
  d <- set$dat |> 
    filter(id %in% prots$id) |> 
    mutate(val = get(what)) |> 
    left_join(set$metadata, by = "sample") |> 
    filter(!(treatment == "Neg")) |> 
    droplevels() |> 
    left_join(prots, by = "id") |> 
    drop_na() |> 
    group_by(id, gene_symbol, group, treatment, time_point) |> 
    summarise(val = mean(val, na.rm = TRUE)) |> 
    ungroup() |> 
    group_by(id) |> 
    mutate(M = mean(val, na.rm = TRUE)) |> 
    mutate(val = (val - M) / log10(2)) |> 
    ungroup() |> 
    mutate(gene_symbol = factor(gene_symbol, levels = genes) |> fct_rev()) |> 
    ungroup() |> 
    add_axis_groups("group", "treatment")
  

  if (is.null(max_z)) 
    max_z <- max(abs(d$val))
  
  
  d_time <- d |> 
    select(x, time_point) |> 
    distinct()
  d_treat <- d |> 
    select(x, treatment) |> 
    distinct() |> 
    group_by(treatment) |> 
    summarise(x = mean(x))
  
  g1 <- ggplot() +
    theme_void() +
    geom_text(data = d_time, aes(x = x, y = 1, label = time_point)) +
    geom_text(data = d_treat, aes(x = x, y = 2, label = treatment)) +
    scale_x_continuous(expand = c(0, 0), limits = c(0.5, 9.5)) +
    theme(plot.margin = margin(5, 0, 5, 0))
    
  g2 <- d |> 
    ggplot(aes(x = x, y = gene_symbol, fill = val)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      text = element_text(size = text_size),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    geom_tile() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = NULL, y = NULL, fill = expression(log[2]~FC)) +
    scale_fill_distiller(type = "div", palette = "RdBu", limits = c(-max_z, max_z), expand = c(0, 0))
  
  #cowplot::plot_grid(g1, g2, ncol = 1, rel_heights = c(1, 10), align = "v", axis = "lr")
  g1 / g2 +
    plot_layout(heights = c(1, 10))
    
}

mn_plot_volcano <- function(da, ctr, pgr = NULL, fdr_limit = 0.01, logfc_limit = 0, with_names = TRUE, palette = okabe_ito_palette) {
  d <- da |> 
    filter(contrast == ctr) |> 
    mutate(
      name = str_remove(gene_symbols, ";.*$")
    )
  
  if(is.null(pgr)) {
    d <- d |> 
      mutate(
        sig = FDR < fdr_limit & abs(logFC) > logfc_limit,
        sel = sig,
        group = sig
      )
  } else {
    d <- d |> 
      left_join(pgr, by = "id") |> 
      mutate(sel = !is.na(group))
  }
  
  d_nsel <- d |> filter(!sel)
  d_sel <- d |> filter(sel)
  
  g <- ggplot(d, aes(x = logFC, y = -log10(PValue), label = name)) +
    th +
    geom_point(data = d_nsel, colour = "grey80", size = 0.8) +
    geom_point(data = d_sel, aes(colour = group), size = 1.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    #geom_vline(xintercept = c(-logfc_limit, logfc_limit), linetype = "dashed", linewidth = 0.2, colour = "darkgreen") +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P))
  
  if(is.null(pgr)) {
    g <- g +
      theme(legend.position = "none") +
      scale_colour_manual(values = "black")
  } else {
    g <- g +
      scale_colour_manual(values = palette)
  }
  
  if(with_names) {
    g <- g + geom_text_repel(data = d_sel, size = 2.5, max.overlaps = 100) 
  }
  
  sigfilt_01 <- d |> filter(FDR < 0.01)
  sigfilt_05 <- d |> filter(FDR < 0.05)
  if(nrow(sigfilt_01) > 0) {
    siglim <- -log10(sigfilt_01$PValue) |> min()
    g <- g + geom_hline(yintercept = siglim, linetype = "dashed", linewidth = 0.2, colour = "darkgreen")
  }
  if(nrow(sigfilt_05) > 0) {
    siglim <- -log10(sigfilt_05$PValue) |> min()
    g <- g + geom_hline(yintercept = siglim, linetype = "dashed", linewidth = 0.2, colour = "darkgreen")
  }
  g
}



mn_plot_ma <- function(da, ctr, pgr = NULL, fdr_limit = 0.01, logfc_limit = 0, with_names = TRUE, palette = okabe_ito_palette) {
  d <- da |> 
    filter(contrast == ctr) |> 
    mutate(
      name = str_remove(gene_symbols, ";.*$")
    )
  
  if(is.null(pgr)) {
    d <- d |> 
      mutate(
        sig = FDR < fdr_limit & abs(logFC) > logfc_limit,
        sel = sig,
        group = sig
      )
  } else {
    d <- d |> 
      left_join(pgr, by = "id") |> 
      mutate(sel = !is.na(group))
  }
  
  d_nsel <- d |> filter(!sel)
  d_sel <- d |> filter(sel)
  
  g <- ggplot(d, aes(x = AveExpr, y = logFC, label = name)) +
    th +
    geom_point(data = d_nsel, colour = "grey80", size = 0.8) +
    geom_point(data = d_sel, aes(colour = group), size = 1.3) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(y = expression(log[2]~FC), x = expression(log[10]~Intensity))
  
  if(is.null(pgr)) {
    g <- g +
      theme(legend.position = "none") +
      scale_colour_manual(values = "black")
  } else {
    g <- g +
      scale_colour_manual(values = palette)
  }
  
  if(with_names) {
    g <- g + geom_text_repel(data = d_sel, size = 2.5, max.overlaps = 100) 
  }
  
  g
}


mn_plot_volma <- function(da, ctr, name, pgr = NULL, fdr_limit = 0.01, logfc_limit = 0, with_names = TRUE, palette = okabe_ito_palette) {
  mn_plot_volcano(da, ctr, pgr, fdr_limit, logfc_limit, with_names, palette) |> 
    gp(str_glue("volcano_{name}"), 6, 6)
  mn_plot_ma(da, ctr, pgr, fdr_limit, logfc_limit, with_names, palette) |> 
    gp(str_glue("ma_{name}"), 6, 6)
}


read_proteins_for_heatmaps <- function(prot_file) {
  sheets <- excel_sheets(prot_file)
  map(sheets, function(sheet) {
    read_excel(prot_file, sheet = sheet) |> 
      add_column(figure = sheet)
  }) |> 
    list_rbind()
}


get_heatmap_data <- function(set, genes = NULL, where = "", what = "abu_limma", mean_condition = TRUE) {
  prots <- set$id_prot_gene |> 
    select(id, gene_symbol = gene_symbols) |> 
    distinct()
  
  if(!is.null(genes)) {
    prots <- prots |> 
      filter(gene_symbol %in% genes) 
  }
  
  if(nrow(prots) != length(genes)) {
    print(paste("Missing proteins in", where, setdiff(genes, prots$gene_symbol) |> paste(collapse = ", ")))
  }
    
  
  d <- set$dat |> 
    filter(id %in% prots$id) |> 
    mutate(val = get(what)) |> 
    left_join(set$metadata, by = "sample") |> 
    filter(!(treatment == "Neg")) |> 
    droplevels() |> 
    left_join(prots, by = "id") |> 
    drop_na() |> 
    unite(grp, c(treatment, time_point, replicate), remove = FALSE)
  
  if(mean_condition) {
    d <- d |> 
      group_by(id, gene_symbol, group, treatment, time_point) |> 
      summarise(val = mean(val, na.rm = TRUE)) |> 
      ungroup() |> 
      mutate(grp = group)
  }
  
  d |> 
    group_by(id) |> 
    mutate(M = mean(val, na.rm = TRUE)) |> 
    mutate(val = (val - M) / log10(2)) |> 
    ungroup() |> 
    mutate(gene_symbol = factor(gene_symbol, levels = genes) |> fct_rev()) |> 
    pivot_wider(id_cols = id, names_from = grp, values_from = val) |> 
    left_join(prots, by = "id") |> 
    select(-id) |> 
    relocate(gene_symbol, .before = 1)
}


make_proteins_for_heatmaps <- function(set, hprots, suffix) {
  pth <- str_glue("for_manuscript/proteins_for_heatmaps_{suffix}")
  if(!dir.exists(pth)) dir.create(pth)
  
  figs <- hprots$figure |> unique()
  for(f in figs) {
    fname = file.path(pth, paste0(f, ".csv"))
    genes <- hprots |> 
      filter(figure == f) |> 
      pull(name) |> 
      unique()
    get_heatmap_data(set, genes, f) |> 
      mutate(across(where(is.numeric), ~signif(.x, digits = 4))) |> 
      write_csv(fname)
  }
}

table_heatmap <- function(set, what = "abu_limma") {
  meta <- set$metadata |> 
    filter(!(treatment == "Neg")) |> 
    droplevels() 
  
  dat <- set$dat |> 
    mutate(val = get(what)) |> 
    inner_join(meta, by = "sample")
  
  make_heat_ <- function(d) {
    d |> 
      group_by(id) |> 
      mutate(M = mean(val, na.rm = TRUE)) |> 
      mutate(val = (val - M) / log10(2)) |> 
      ungroup() |> 
      pivot_wider(id_cols = id, names_from = group, values_from = val)
  }
  
  d_rep <- dat |> 
    select(-group) |> 
    unite(group, c(treatment, time_point, replicate), remove = FALSE) |> 
    arrange(treatment, time_point, replicate) |> 
    make_heat_()
    
  d_mean <- dat |> 
    group_by(id, group, treatment, time_point) |> 
    summarise(val = mean(val, na.rm = TRUE)) |> 
    ungroup() |> 
    arrange(treatment, time_point) |> 
    make_heat_()
  
  set$info |>
    select(id, proteins = protein, gene_symbols) |> 
    left_join(d_rep, by = "id") |> 
    left_join(d_mean, by = "id")
}


mn_fig_3ab <- function(da) {
  i2g <- set$info |> 
    select(id, gene_symbols, protein)
  
  ip <- i2g |> select(id, protein) |> separate_longer_delim(protein, delim = ";")
  ig <- i2g |> select(id, gene_symbols) |> separate_longer_delim(gene_symbols, delim = ";")
  
  prots <- read_excel("for_manuscript/Figure 3A_3B_DNA repair proteins.xlsx") |> 
    set_names("group", "gene_symbols", "protein") |> 
    drop_na()
  
  
    drop_na() |> 
    left_join(i2g, by = "gene_symbols") |> 
    drop_na()

  mn_plot_volcano(da, "treatmentDRB", pgr, with_names = FALSE) |> gp("fig_3a", 6, 6)
  mn_plot_volcano(da, "treatmentTPL", pgr, with_names = FALSE) |> gp("fig_3b", 6, 6)
}

mn_plot_significant_ <- function(da, fdr_limit = 0.01) {
  n_tot <- da$id |> unique() |> length()
  ids_drb <- da |> filter(contrast == "treatmentDRB" & FDR < fdr_limit) |> pull(id)
  ids_tpl <- da |> filter(contrast == "treatmentTPL" & FDR < fdr_limit) |> pull(id)
  tribble(
    ~what, ~count,
    "DRB", length(ids_drb),
    "TPL", length(ids_tpl),
    "DRB+TPL", length(intersect(ids_drb, ids_tpl))
  ) |> 
    mutate(perc = 100 * count / n_tot) |> 
    ggplot(aes(x = what, y = perc, label = count)) +
    th +
    geom_col() +
    geom_text(vjust = -1) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Percentage")
}


get_interactor_genes <- function(prot_file) {
  sheets <- excel_sheets(prot_file)
  map(sheets, function(sheet) {
    read_excel(prot_file, sheet = sheet) |> 
      janitor::clean_names() |> 
      select(gene_symbol = protein_names) |> 
      add_column(interactor = sheet, group = "RNAPII interactor", group_name = "inter")
  }) |> 
    list_rbind()
}

get_short_lived_genes <- function(prot_file) {
  read_excel(prot_file) |> 
    janitor::clean_names() |> 
    select(protein = uniprot_id, gene_symbol) |> 
    add_column(group = "Short lived", group_name = "short")
}

get_gene_groups <- function() {
  bind_rows(
    get_interactor_genes("for_manuscript/RNAPII interactor_fromFigS6.xlsx"),
    get_short_lived_genes("for_manuscript/Figure 1J_Short half life list_RPE1 cells.xlsx")
  )
}



full_gene_sel <- function(da, sel = NULL) {
  if(!is.null(sel)) {
    gsel <- da |> 
      select(id, gene_symbols) |> 
      separate_longer_delim(gene_symbols, delim = ";") |> 
      distinct() |> 
      filter(gene_symbols %in% sel) |> 
      pull(id)
    
    da <- da |>
      filter(id %in% gsel)
  }
  da
}

count_significant <- function(da, sel = NULL, sel_are_ids, fdr_limit = 0.01) {
  
  if(sel_are_ids) {
    da <- da |>
      filter(id %in% sel)
  } else {
    da <- full_gene_sel(da, sel)
  }
  
  ds <- da |> 
    filter(FDR < fdr_limit & contrast %in% c("treatmentDRB", "treatmentTPL"))
  
  n_tot <- da$id |> unique() |> length()
  n_sig <- ds$id |> unique() |> length()
  
  genes_drb <- ds |> filter(contrast == "treatmentDRB") |> pull(gene_symbols)
  genes_tpl <- ds |> filter(contrast == "treatmentTPL") |> pull(gene_symbols)
  tribble(
    ~what, ~count,
    "DRB", length(genes_drb),
    "TPL", length(genes_tpl),
    "DRB+TPL", length(intersect(genes_drb, genes_tpl)),
    "Not affected", n_tot - n_sig,
    "Total", n_tot
  ) |> 
    mutate(perc = 100 * count / n_tot) |> 
    add_column(total = n_tot) |> 
    nest(data = c(count, total)) |> 
    mutate(
      fit = map(data, ~prop.test(.x$count, .x$total)),
      tidied = map(fit, broom::tidy)
    ) |> 
    unnest(c(data, tidied)) |> 
    select(-c(fit, method, alternative)) |> 
    mutate(conf.low = 100 * conf.low, conf.high = 100 * conf.high)
}

mn_find_sig_sel <- function(da, sel) {
  full_gene_sel(da, sel) |> 
    filter(sig & str_detect(contrast, "treat")) |> 
    select(gene_symbols, contrast)
}


mn_plot_significant <- function(da, title, suffix, sel = NULL, sel_are_ids = FALSE) {
  
  sel <- unique(sel)
  
  pl <- function(fdr_limit) {
    cs <- count_significant(da, sel, sel_are_ids, fdr_limit = fdr_limit)
    tot <- cs |> filter(what == "Total") |> pull(count)
    
    g <- cs |> 
      filter(!(what %in% c("Not affected", "Total"))) |> 
      ggplot(aes(x = what, y = perc, label = count, fill = what)) +
      th +
      theme(legend.position = "none") +
      geom_col(colour = "grey30") +
      geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
      geom_text(vjust = -1, hjust = 1.3) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
      scale_fill_manual(values = okabe_ito_palette) +
      labs(x = NULL, y = "Percentage", fill = "Selection", title = str_glue("{title} ({tot})"))
    gp(g, str_glue("significant_drb_tpl_{suffix}_{fdr_limit}"), 3, 3)
  }
  
  pl(0.01)
  pl(0.05)
}

mn_plot_significant_groups <- function(da, gene_groups, name) {
  nm <- str_replace(name, "_", " ")
  mn_plot_significant(da, str_glue("All proteins {nm}"), str_glue("all_{name}"))

  gene_groups |> 
    group_split(group, group_name) |> 
    map(function(w) {
      # Bottom 20% are ids, the rest are gene_symbols
      if(all(is.na(w$id))) {
        sel <- w$gene_symbol
        sel_are_ids <- FALSE
      } else {
        sel <- w$id
        sel_are_ids <- TRUE
      }
      mn_plot_significant(da, paste(w$group[1], nm), paste0(w$group_name[1], "_", name), sel, sel_are_ids)
    })
}

dna_repair_pathways <- c(
  "GO:0006281",
  "GO:0006282",
  "GO:0006298",
  "GO:0006284",
  "GO:0006289"
)

dna_repair <- "GO:0006281"

mn_volcano_pathways <- function(da, ctr, terms, pathways) {
  trm <- terms$go$mapping |> 
    filter(term_id %in% pathways) |> 
    select(gene_symbol, term_id) |> 
    distinct() |> 
    left_join(terms$go$terms, by = "term_id")
  
  pgr <- da |> 
    get_ids() |> 
    right_join(trm, by = "gene_symbol") |> 
    select(id, group = term_name) |> 
    distinct()
  
  mn_plot_volcano(da, ctr, pgr, with_names = FALSE, palette = "black") + theme(legend.position = "none")
}

# genes is a tibble with gene_symbol and group
mn_volcano_groups<- function(da, ctr, genes, with_names = TRUE) {
  pgr <- da |> 
    get_ids() |> 
    right_join(select(genes, gene_symbol, group), by = "gene_symbol", relationship = "many-to-many") |> 
    select(id, group) |> 
    drop_na() |> 
    distinct()
  
  d <- da |> 
    filter(contrast == ctr) |> 
    mutate(name = str_remove(gene_symbols, ";.*$")) |> 
    left_join(pgr, by = "id") |> 
    mutate(sel = !is.na(group))

  d_nsel <- d |> filter(!sel)
  d_sel <- d |> filter(sel)
  
  g <- ggplot(d, aes(x = logFC, y = -log10(PValue), label = name)) +
    th +
    geom_point(data = d_nsel, colour = "grey80", size = 0.8) +
    geom_point(data = d_sel, colour = "black", size = 1.3) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    facet_wrap(~group)
  
  if(with_names) {
    g <- g + geom_text_repel(data = d_sel, size = 2.5, max.overlaps = 100) 
  }
  g  
}


get_tfs <- function(prot_file) {
  read_excel(prot_file) |> 
    janitor::clean_names() |> 
    select(protein = uni_prot_id, gene_symbol = hgnc_approved_gene_symbol)
}

mn_volcano_tfs <- function(da, ctr) {
  tfs <- get_tfs("for_manuscript/TF_List_LambertAnddb.xlsx") |> 
    add_column(group = "transcription factor")
  
  pgr <- da |> 
    get_ids() |> 
    right_join(tfs, by = "gene_symbol") |> 
    select(id, group) |> 
    distinct()
  
  mn_plot_volcano(da, ctr, pgr, with_names = FALSE, palette = "black") + theme(legend.position = "none")
}


mn_plot_transfac_counts <- function() {
  d <- read_excel("for_manuscript/transcription_factor_counts.xlsx") |> 
    mutate(total = sum(count)) |> 
    mutate(perc = 100 * count / total) |> 
    nest(data = c(count, total)) |> 
    mutate(
      fit = map(data, ~prop.test(.x$count, .x$total)),
      tidied = map(fit, broom::tidy)
    ) |> 
    unnest(c(data, tidied)) |> 
    select(-c(fit, method, alternative)) |> 
    mutate(conf.low = 100 * conf.low, conf.high = 100 * conf.high) |> 
    mutate(what = as_factor(what) |> fct_rev())
  #tot <- sum(d$count)
  d |> 
    ggplot(aes(x = what, y = perc, label = count, fill = what)) +
    th +
    theme(
      legend.position = "none"
      #axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    geom_col(colour = "grey30") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1) +
    geom_text(vjust = -0.8, hjust = -0.2) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 15)) +
    scale_fill_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Percentage", fill = "Selection") +
    coord_flip()
}


mean_tp <- function(set, what = "abu_limma") {
  set$dat |> 
    mutate(val = get(what)) |> 
    left_join(set$metadata, by = "sample") |> 
    filter(!(treatment == "Neg")) |> 
    droplevels() |> 
    group_by(id, group) |> 
    summarise(val = mean(val, na.rm = TRUE)) |> 
    ungroup() |> 
    pivot_wider(id_cols = id, names_from = grp, values_from = val) |> 
    column_to_rownames("i")
    as.matrix()
}

mn_ip_input_correlation <- function(set_ip, set_input, what = "abu_limma") {
  set <- merge_sets(set_ip, set_input)
  meta <- set$metadata |> 
    filter(!(treatment == "Neg")) |> 
    unite(grp, c(protocol, treatment, time_point), remove = FALSE) |> 
    droplevels()

  M <- set$dat |> 
    mutate(val = get(what)) |> 
    right_join(meta, by = "sample") |> 
    group_by(id, grp) |> 
    summarise(val = mean(val, na.rm = TRUE)) |> 
    ungroup() |> 
    pivot_wider(id_cols = id, names_from = grp, values_from = val) |> 
    column_to_rownames("id") |> 
    as.matrix()
  
  m <- meta |> 
    select(protocol, grp) |> 
    distinct() |>
    mutate(grp = as.character(grp))
  
  M_ip <- M[, m |> filter(protocol == "IP") |> pull(grp)]
  M_inp <- M[, m |> filter(protocol == "Input") |> pull(grp)]
  
  cor(M_ip, M_inp, use = "pairwise.complete.obs") |> 
    as.data.frame() |> 
    rownames_to_column("p_ip") |> 
    pivot_longer(-p_ip, names_to = "p_inp") |> 
    ggplot(aes(x = p_ip, y = p_inp, fill = value)) +
    th +
    geom_tile() +
    scale_fill_viridis_c() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(x = NULL, y = NULL)
}

read_protein_groups_id <- function(efile, da) {
  genes <- read_excel(efile) |> 
    select(group, gene_symbol) |> 
    distinct()
  
  gene2id <- da |> 
    get_ids() |> 
    right_join(genes, by = "gene_symbol", relationship = "many-to-many")
  
  print("Missing:")
  gene2id |> filter(is.na(id)) |> select(group, gene_symbol) |> print(n = 1000)
  
  gene2id |> 
    drop_na()
}


gse_boot <- function(da, ctr, treat, expmt, protcl, genes, n_boot = 1000, seed = 42) {
  ds <- da |> 
    filter(contrast == ctr)
  
  all_ids <- unique(ds$id)
  mean_fdr <- function(sel) {
    mean(ds[ds$id %in% sel, ]$logFC)
  }
  
  do_boot <- function(data) {
    n_obs <- nrow(data)
    m_obs <- mean_fdr(data$id)
    boot <- map(1:n_boot, function(i) {
      sel <- sample(all_ids, n_obs)
      mean_fdr(sel)
    }) |> 
      unlist()
    tibble(
      n_obs = n_obs,
      m_obs = m_obs,
      n_boot = n_boot,
      n_less = length(which(boot < m_obs)),
      p_boot =  n_less / n_boot,
      contrast = ctr,
      experiment = expmt,
      protocol = protcl,
      treatment = treat,
    )
  }
  
  do_boot_p <- function(data) {
    n_obs <- nrow(data)
    m_obs <- mean_fdr(data$id)
    n_chunk <- 5
    n_subboot <- n_boot / n_chunk
    boot <- parallel::mclapply(1:n_chunk, function(i) {
      map(1:n_subboot, function(k) {
        sel <- sample(all_ids, n_obs)
        mean_fdr(sel)
      })
    }, mc.cores = 6) |> 
      unlist()
    tibble(
      n_obs = n_obs,
      m_obs = m_obs,
      n_boot = n_boot,
      n_less = length(which(boot < m_obs)),
      p_boot =  n_less / n_boot,
      contrast = ctr,
      experiment = expmt,
      protocol = protcl,
      treatment = treat,
    )
  }
  
  set.seed(seed)
  genes |> 
    select(group, id, gene_symbol) |> 
    nest(data = c(id, gene_symbol)) |> 
    mutate(res = map(data, ~do_boot_p(.x))) |> 
    unnest(res) |> 
    select(group, experiment, protocol, treatment, size = n_obs, mean_logFC = m_obs, n_boot, n_less, p_value = p_boot)
}


mn_plot_group_boot <- function(mdl) {
  
  m <- mdl |>
    unite(expm, c(experiment, protocol), sep = " ") |> 
    mutate(expm = fct_relevel(expm, c("E1 IP", "E2 IP"))) |> 
    group_by(expm) |> 
    mutate(
      fdr = p.adjust(p_value, method = "BH"),
      sig = fdr < 0.05
    )
  ms <- m |> 
    filter(sig) |> 
    mutate(ast = (abs(mean_logFC) + 0.03) * sign(mean_logFC))
  
  m |> 
    ggplot(aes(x = group, y = mean_logFC)) +
    th +
    theme(legend.position = "right") +
    geom_col(aes(fill = log10(fdr))) +
    geom_hline(yintercept = 0, colour = "grey50") +
    geom_point(data = ms, aes(y = ast), shape = 8, size = 0.8) +
    scale_fill_viridis_c(option = "cividis") +
    facet_grid(expm ~ treatment) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05)), breaks = c(-0.4, -0.2, 0, 0.2)) +
    coord_flip() +
    labs(y = expression(mean~log[2]~FC), x = NULL, fill = expression(log[10]~P))
}

mn_venn_input_ip <- function(da_input, da_ip, id_mapping) {
  id_map <- id_mapping |> 
    drop_na() |> 
    filter(n.x == 1 & n.y == 1) |> 
    select(-c(n.x, n.y))
  
  # significant in any treatment
  get_sig_ <- function(da) {
    da |> 
      mutate(sig = FDR < 0.05) |> 
      filter(str_detect(contrast, "treatment")) |> 
      select(id, contrast, sig) |> 
      pivot_wider(id_cols = id, names_from = contrast, values_from = sig) |> 
      mutate(sig = treatmentDRB | treatmentTPL) |> 
      select(id, sig)
  }
  
  ds_ip <- get_sig_(da_ip)
  ds_input <- get_sig_(da_input)

  ds_input |> 
    inner_join(id_map, by = c("id" = "id.x")) |> 
    left_join(ds_ip, by = c("id.y" = "id")) |> 
    drop_na() |> 
    select(Input = sig.x, IP = sig.y, id.input = id, id.ip = id.y)
}


mn_venn_groups <- function(da, gene_groups) {
  ds <- da |> 
    filter(str_detect(contrast, "treatment")) |> 
    mutate(sig = FDR < 0.05) |> 
    select(id, sig, gene_symbols) |> 
    distinct()
  
  genes_grp <- gene_groups |> 
    group_split(group, group_name) |> 
    map(function(w) {
      # Bottom 20% are ids, the rest are gene_symbols
      if(!all(is.na(w$id))) {
        dg <- ds |> filter(id %in% w$id)
      } else {
        dg <- full_gene_sel(ds, w$gene_symbol)
      }
      dg |> 
        select(id) |> 
        distinct() |> 
        add_column(group = w$group[1])
    }) |> 
    list_rbind()
  genes_sig <- ds |> 
    filter(sig) |> 
    select(id) |> 
    distinct() |> 
    add_column(group = "Significant")
  
  bind_rows(genes_sig, genes_grp) |> 
    mutate(val = TRUE) |> 
    pivot_wider(id_cols = id, names_from = group, values_from = val, values_fill = FALSE) |> 
    select(-id)
}


groups_for_gsea <- function(gene_grp) {
  groups <- gene_grp$group |> unique()
  term2feature <- map(groups, function(grp) {
    gene_grp |> 
      filter(group == grp) |> 
      pull(id)
  }) |> 
    set_names(groups)
  
  term2name <- list(groups = groups)
  
  list(
    group = list(
      term2feature = term2feature,
      term2name = term2name
    )
  )
}


scientific_10 <- function(x) {
  parse(
    text = scales::scientific_format()(x) |>
      str_replace("e\\+*", " %*% 10^") |>
      str_replace("^1 \\%\\*\\% 10\\^", "10^")
  )
}


plot_fig_1d <- function(data_file) {
  pl_ <- function(sheet) {
    d <- read_excel(data_file, sheet = sheet) |> 
      pivot_longer(everything()) |> 
      drop_na() |> 
      mutate(name = fct_relevel(name, c("no EdU", "DMSO")))
    md <- d |>
      group_by(name) |> 
      summarise(m = median(value)) |> 
      mutate(xi = as.integer(name))
    dd <- 0.4
    
    tst <- wilcox.test(value ~ name, data = d |> filter(name != "no EdU"))
    print(tst$p.value)
    
    pl <- d |> 
      ggplot(aes(x = name, y = value)) +
      theme_classic() +
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)
      ) +
      geom_quasirandom(size = 0.05) +
      scale_y_log10(labels = scientific_10, expand = c(0, 0), limits = c(10, 1e5)) +
      labs(x = NULL, y = "EdU mean intensity in EdU + cells (a.u.)") +
      geom_segment(data = md, aes(x = xi - dd, xend = xi + dd, y = m, yend = m),
                   colour = "red", alpha = 0.5, linewidth = 1.3)
    gp(pl, sheet, 1.8, 3.2)
  }
  
  pl_("Fig1D_TPL")
  pl_("Fig1D_DRB")
  
  
}