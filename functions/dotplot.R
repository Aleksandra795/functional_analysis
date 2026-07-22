# CERNO dot plot

cerno_dotplot <- function(df,
                          top_n = 25,
                          padj_cutoff = 0.05,
                          use_top_n_if_no_significant = TRUE,
                          id_col = "ID",
                          term_col = "Title",
                          padj_col = "adj.P.Val",
                          auc_col = "AUC",
                          highlight1 = NULL,
                          highlight2 = NULL,
                          highlight_color1 = "#009E73",
                          highlight_color2 = "#0072B2",
                          split_after_words = 7,
                          point_alpha = 0.9,
                          point_size_range = c(2.5, 9)) {

  required_cols <- c(term_col, padj_col, auc_col)
  missing_cols <- setdiff(required_cols, colnames(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  plot_df <- df %>%
    mutate(
      pathway_raw = .data[[term_col]],
      pathway_clean = clean_pathway_name(
        .data[[term_col]],
        split_after_words = split_after_words
      ),
      neglog10_padj = -log10(pmax(.data[[padj_col]], .Machine$double.xmin))
    ) %>%
    filter(
      !is.na(.data[[padj_col]]),
      !is.na(.data[[auc_col]])
    ) %>%
    arrange(.data[[padj_col]])

  significant_df <- plot_df %>%
    filter(.data[[padj_col]] < padj_cutoff)

  if (nrow(significant_df) > 0) {
    plot_df <- significant_df
  } else if (isTRUE(use_top_n_if_no_significant)) {
    plot_df <- plot_df %>%
      slice_head(n = top_n)
  } else {
    stop("No pathways meet padj_cutoff.")
  }

  plot_df <- plot_df %>%
    arrange(.data[[padj_col]]) %>%
    slice_head(n = top_n) %>%
    mutate(
      pathway_clean = factor(
        pathway_clean,
        levels = rev(unique(pathway_clean))
      )
    )

  label_levels <- levels(plot_df$pathway_clean)

  label_map <- setNames(
    highlight_pathway_labels(
      labels = label_levels,
      highlight1 = highlight1,
      highlight2 = highlight2,
      color1 = highlight_color1,
      color2 = highlight_color2
    ),
    label_levels
  )

  ggplot(
    plot_df,
    aes(
      x = neglog10_padj,
      y = pathway_clean,
      size = .data[[auc_col]],
      color = neglog10_padj
    )
  ) +
    geom_point(alpha = point_alpha) +
    scale_y_discrete(labels = label_map) +
    scale_size_continuous(
      range = point_size_range,
      name = "AUC"
    ) +
    labs(
      x = expression(-log[10]("FDR")),
      y = NULL,
      color = expression(-log[10]("FDR"))
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.y = ggtext::element_markdown(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )
}
