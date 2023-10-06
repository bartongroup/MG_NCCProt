okabe_ito_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "grey80", "grey30", "black")

gs <- function(gg, name, width, height) {
  ggsave(filename = file.path("fig", paste0(name, ".png")), plot = gg, device = "png",
         width = width, height = height, dpi = 300)
}



plot_detection <- function(set) {
  cumcurve <- function(d, what) {
    d |>
      group_by(get(what)) |> 
      tally() |>
      arrange(desc(n)) |> 
      mutate(x = row_number())
  }
  plotcurve <- function(d) {
    ggplot(d, aes(x = x, y = n)) +
      theme_bw() +
      geom_step(direction = "vh") +
      scale_x_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA)) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.03)), limits = c(0, NA))
  }
  
  dp <- cumcurve(set$dat, "id")
  ds <- cumcurve(set$dat, "sample")
  
  rm(set)
  
  g1 <- plotcurve(dp) + labs(x = "Proteins", y = "Detected in that many samples")
  g2 <- plotcurve(ds) + labs(x = "Samples", y = "Contains that many proteins")
  
  cowplot::plot_grid(g1, g2, align = "h", labels = c("A", "B"))
}

plot_sample_detection <- function(set) {
  d <- set$dat |> 
    group_by(sample) |> 
    tally()
  rm(set)
  ggplot(d, aes(x = n, y = fct_reorder(sample, n))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_segment(aes(xend = 0, yend = fct_reorder(sample, n)), colour = "grey70") +
    geom_point() +
    labs(x = "Number of detected proteins", y = "Sample")
}


plot_sample_ridges <- function(st, what = "abu_med", name_var = "sample", fill_var = "treatment", scale = 3, bandwidth = 0.1) {
  env <- new.env(parent = globalenv())
  env$fill_var <- fill_var
  env$scale <- scale
  env$what <- what
  env$name_var <- name_var
  env$bandwidth <- bandwidth
  
  levs <- st$dat |>
    rename(gr := !!name_var) |> 
    group_by(gr) |>
    summarise(m = median(get(what))) |>
    arrange(m) |>
    pull(gr) 
  
  env$d <- st$dat |>
    left_join(st$metadata, name_var) |>
    mutate(gr = factor(get(name_var), levels = levs)) |> 
    mutate(val = get(what), fil = get(fill_var)) |> 
    select(gr, val, fil)
  
  with(env, {
    d |> 
      ggplot(aes(x = val, y = gr, fill = fil)) +
      theme_bw() +
      theme(panel.grid = element_blank()) +
      geom_density_ridges(scale = scale, bandwidth = bandwidth) +
      #geom_vline(xintercept = 0) +
      scale_fill_manual(values = okabe_ito_palette) +
      labs(x = what, y = name_var, fill = fill_var)    
  })
}


plot_clustering <- function(set, text_size = 10, what = "abu_med", dist.method = "euclidean",
                            clust.method = "complete", colour_var = "treatment") {
  tab <- dat2mat(set$dat, what)
  
  dendr <- t(tab) |> 
    dist(method = dist.method) |> 
    hclust(method = clust.method) |>
    dendsort::dendsort() |> 
    ggdendro::dendro_data()
  
  seg <- ggdendro::segment(dendr)
  meta <- set$metadata |>
    mutate(
      colvar = get(colour_var),
      sample = as.character(sample)
    )
  labs <- left_join(dendr$labels |> mutate(label = as.character(label)), meta, by = c("label" = "sample")) |> 
    mutate(colour = okabe_ito_palette[as_factor(colvar)])
  theme.d <- ggplot2::theme(
    panel.grid.major = ggplot2::element_blank(),
    panel.grid.minor = ggplot2::element_blank(),
    panel.background = ggplot2::element_blank(),
    axis.text.y = ggplot2::element_text(size = text_size, colour = labs$colour),
    axis.line.y = ggplot2::element_blank(),
    axis.line.x = ggplot2::element_line(linewidth = 0.5),
    axis.ticks.y = ggplot2::element_blank()
  )
  rm(set, tab, dendr)
  ggplot() +
    theme.d +
    coord_flip() +
    geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend)) +
    scale_x_continuous(breaks = seq_along(labs$label), labels = labs$label) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(seg$y) * 1.03)) +
    scale_colour_manual(values = okabe_ito_palette) +
    labs(x = NULL, y = "Distance")
}


plot_distance_matrix <- function(set, what = "abu_med", text_size = 10) {
  tab <- dat2mat(set$dat, what)
  
  d <- cor(tab, use = "complete.obs") |> 
    as_tibble(rownames = "sample") |>
    pivot_longer(-sample) |> 
    mutate(sample = factor(sample, levels = set$metadata$sample)) |> 
    mutate(name = factor(name, levels = set$metadata$sample))
  
  rm(set, tab)
  ggplot(d, aes(x = sample, y = name)) +
    geom_tile(aes(fill = value)) +
    scale_fill_viridis_c(option = "cividis") +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = text_size),
      axis.text.y = element_text(size = text_size),
      legend.position = "top"
    ) +
    labs(x = NULL, y = NULL, fill = "Correlation")
}


plot_xy <- function(dat, colour_var, shape_var, point_size = 1) {
  dat |> 
    rename(colvar = !!colour_var, shapevar = !!shape_var) |> 
    ggplot(aes(x = x, y = y, colour = colvar, shape = shapevar)) +
    theme_bw() +
    geom_point(size = point_size) +
    scale_color_manual(values = okabe_ito_palette, name = colour_var) +
    scale_shape_manual(values = c(15:18, 0, 1, 2, 5, 6), na.value = 4, name = shape_var) +
    labs(x = NULL, y = NULL) +
    theme(
      panel.grid = element_blank()
    ) 
  #geom_text_repel(aes(label = sample), size = 3, colour = "black")
}

pca2xy <- function(pc, meta) {
  pc$x |> 
    as_tibble(rownames = "sample") |> 
    select(x = PC1, y = PC2, sample) |> 
    left_join(meta, by = "sample")
}


plot_pca <- function(set, point_size = 2, what = "abu_med", colour_var = "treatment", shape_var = "time_point") {
  tab <- dat2mat(set$dat, what)
  
  # remove rows with zero variance
  tab <- tab[apply(tab, 1, function(v) sum(is.na(v)) == 0), ]
  tab <- tab[apply(tab, 1, sd) > 0, ]
  
  pca <- prcomp(t(tab), scale. = TRUE, center = TRUE)
  var.perc <- 100 * (pca$sdev)^2 / sum((pca$sdev)^2)
  pca1 <- sprintf("PCA1 (%5.1f%%)", var.perc[1])
  pca2 <- sprintf("PCA2 (%5.1f%%)", var.perc[2])
  pca2xy(pca, set$metadata) |> 
    plot_xy(colour_var, shape_var, point_size) +
    labs(x = pca1, y = pca2)
}

plot_volma <- function(res, group, point_size, point_alpha) {
  r <- res |>
    mutate(
      group = get(group)
    ) |> 
    select(x, y, sel, group)
  r_sel <- r |> filter(sel)
  r_nsel <- r |> filter(!sel)
  
  rm(res, r)  # Minimise environment for serialisation
  
  g <- ggplot() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_point(data = r_nsel, aes(x = x, y = y), colour = "grey70",
               size = point_size, alpha = point_alpha) +
    facet_wrap( ~ group) 
  
  if (nrow(r_sel) > 0) {
    g <- g + geom_point(data = r_sel, aes(x = x, y = y), colour = "black",
                        size = point_size, alpha = point_alpha)
  }
  g
}

plot_ma <- function(res, a = "AveExpr", fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                    fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res |> 
    mutate(
      x = get(a),
      y = get(fc),
      sel = get(fdr) < fdr_limit
    ) |> 
    plot_volma(group, point_size, point_alpha) +
    geom_hline(yintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[10]~Intensity), y = expression(log[2]~FC))
}

plot_volcano <- function(res, fc = "logFC", p = "PValue", fdr = "FDR", group = "contrast",
                         fdr_limit = 0.05, point_size = 0.5, point_alpha = 0.5) {
  res |> 
    mutate(
      x = get(fc),
      y = -log10(get(p)),
      sel = get(fdr) < fdr_limit
    ) |> 
    plot_volma(group, point_size, point_alpha) +
    geom_vline(xintercept = 0, linewidth = 0.1, alpha = 0.5) +
    labs(x = expression(log[2]~FC), y = expression(-log[10]~P)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03)))
}

plot_pdist <- function(res, p = "PValue", group = "contrast", n_bins = 50) {
  brks <- seq(0, 1, length.out = n_bins)
  r <- res |> 
    mutate(p = get(p), grp = get(group)) |> 
    select(p, grp)
  rm(res)
  ggplot(r, aes(x = p, y = after_stat(density))) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_histogram(breaks = brks) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.03))) +
    facet_wrap(~grp)
}


plot_up_down <- function(res, fc = "logFC", fdr = "FDR", group = "contrast", fdr_limit = 0.05) {
  r <- res |>
    mutate(
      direction = if_else(get(fc) > 0, "up", "down"),
      sig = get(fdr) < fdr_limit,
      group = get(group)
    ) |> 
    filter(sig) |> 
    group_by(group, direction, .drop = FALSE) |> 
    tally() |>
    mutate(x = if_else(direction == "down", -n, n)) |> 
    mutate(adj = -1.3 * sign(x) / 2 + 0.5)
  rm(res)
  ggplot(r, aes(x = x, y = group, fill = direction)) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.position = "none"
    ) +
    geom_vline(xintercept = 0, colour = "grey50") +
    geom_col() +
    scale_fill_manual(values = okabe_ito_palette) +
    geom_text(aes(label = n, hjust = adj)) +
    scale_x_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    labs(x = "Count", y = "Coefficient")
}

