# CERNO Manhattan-style plot for GO ontologies

cerno_manhattan_go <- function(df,
                               ontology_col = "ONTOLOGY",
                               term_col = "Title",
                               padj_col = "adj.P.Val",
                               size_col = "AUC",
                               padj_cutoff = 0.05,
                               label_top_n = 12,
                               label_padj_cutoff = NULL,
                               max_terms_per_ontology = NULL,
                               ontology_order = c("BP", "CC", "MF"),
                               split_after_words = 8,
                               point_alpha_significant = 0.95,
                               point_alpha_nonsignificant = 0.30,
                               point_size_range = c(1.8, 7),
                               label_size = 3,
                               jitter_width = 0.28,
                               seed = 123) {

  required_cols <- c(ontology_col, term_col, padj_col, size_col)
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  plot_df <- df %>%
    filter(
      !is.na(.data[[ontology_col]]),
      !is.na(.data[[term_col]]),
      !is.na(.data[[padj_col]]),
      !is.na(.data[[size_col]])
    ) %>%
    mutate(
      ONTOLOGY_plot = as.character(.data[[ontology_col]]),
      ONTOLOGY_plot = factor(
        ONTOLOGY_plot,
        levels = ontology_order[ontology_order %in% unique(ONTOLOGY_plot)]
      ),
      term_clean = clean_pathway_name(
        .data[[term_col]],
        split_after_words = split_after_words
      ),
      neglog10_padj = -log10(pmax(.data[[padj_col]], .Machine$double.xmin)),
      significant = .data[[padj_col]] < padj_cutoff
    ) %>%
    filter(!is.na(ONTOLOGY_plot))

  if (!is.null(max_terms_per_ontology)) {
    plot_df <- plot_df %>%
      group_by(ONTOLOGY_plot) %>%
      arrange(.data[[padj_col]], .by_group = TRUE) %>%
      slice_head(n = max_terms_per_ontology) %>%
      ungroup()
  }

  set.seed(seed)

  plot_df <- plot_df %>%
    mutate(
      x_base = as.numeric(ONTOLOGY_plot),
      x_jitter = x_base + runif(
        n = n(),
        min = -jitter_width,
        max = jitter_width
      )
    )

  if (!is.null(label_padj_cutoff)) {
    label_df <- plot_df %>%
      filter(.data[[padj_col]] < label_padj_cutoff)
  } else {
    label_df <- plot_df %>%
      arrange(.data[[padj_col]]) %>%
      slice_head(n = label_top_n)
  }

  ggplot(
    plot_df,
    aes(
      x = x_jitter,
      y = neglog10_padj
    )
  ) +
    geom_point(
      aes(
        size = .data[[size_col]],
        color = ONTOLOGY_plot,
        alpha = significant
      )
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      linewidth = 0.4
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(
        label = term_clean,
        color = ONTOLOGY_plot
      ),
      size = label_size,
      max.overlaps = Inf,
      min.segment.length = 0,
      box.padding = 0.75,
      point.padding = 0.45,
      force = 3,
      force_pull = 0.15,
      max.time = 5,
      max.iter = 50000,
      segment.alpha = 0.45,
      segment.size = 0.25,
      show.legend = FALSE
    ) +
    scale_x_continuous(
      breaks = seq_along(levels(plot_df$ONTOLOGY_plot)),
      labels = levels(plot_df$ONTOLOGY_plot),
      limits = c(
        1 - jitter_width - 0.15,
        length(levels(plot_df$ONTOLOGY_plot)) + jitter_width + 0.15
      )
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.02, 0.10))
    ) +
    scale_size_continuous(
      range = point_size_range,
      name = size_col
    ) +
    scale_alpha_manual(
      values = c(
        `TRUE` = point_alpha_significant,
        `FALSE` = point_alpha_nonsignificant
      ),
      guide = "none"
    ) +
    labs(
      x = NULL,
      y = expression(-log[10]("FDR")),
      color = "GO ontology"
    ) +
    coord_cartesian(clip = "off") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(face = "bold"),
      legend.position = "right",
      plot.margin = margin(10, 40, 10, 10)
    )
}
