# Hierarchical tree plot based on pathway gene overlap

select_cerno_terms_for_tree <- function(df,
                                        ontology = "BP",
                                        ontology_col = "ONTOLOGY",
                                        padj_col = "adj.P.Val",
                                        top_n = 30,
                                        padj_cutoff = 0.05,
                                        selection_mode = c(
                                          "significant_top_n",
                                          "significant_all",
                                          "top_n",
                                          "significant_or_top_n"
                                        )) {

  selection_mode <- match.arg(selection_mode)

  x <- df |>
    filter(.data[[ontology_col]] == ontology) |>
    filter(!is.na(.data[[padj_col]])) |>
    arrange(.data[[padj_col]])

  switch(
    selection_mode,

    significant_top_n = x |>
      filter(.data[[padj_col]] < padj_cutoff) |>
      slice_head(n = top_n),

    significant_all = x |>
      filter(.data[[padj_col]] < padj_cutoff),

    top_n = x |>
      slice_head(n = top_n),

    significant_or_top_n = {
      sig <- x |>
        filter(.data[[padj_col]] < padj_cutoff)

      if (nrow(sig) >= 3) {
        sig |> slice_head(n = top_n)
      } else {
        x |> slice_head(n = top_n)
      }
    }
  )
}

cerno_treeplot <- function(res_df,
                           members_df,
                           ontology = "BP",
                           ontology_col = "ONTOLOGY",
                           id_col = "ID",
                           term_col = "Title",
                           padj_col = "adj.P.Val",
                           size_col = "AUC",
                           top_n = 30,
                           padj_cutoff = 0.05,
                           selection_mode = c(
                             "significant_top_n",
                             "significant_all",
                             "top_n",
                             "significant_or_top_n"
                           ),
                           nCluster = 5,
                           cluster_method = "ward.D2",
                           branch_length = "none",
                           split_after_words = 7,
                           cluster_label_words = 3,
                           tip_label_size = 2.2,
                           clade_label_size = 2.8,
                           point_size_range = c(0.8, 2.8),
                           tip_label_offset = 0.12,
                           tip_label_space = 6.5,
                           group_label_offset = 3.0,
                           right_padding = 0.8,
                           tree_scale = 0.35,
                           clade_alpha = 0.30,
                           branch_size = 0.28,
                           hilight_extend = 0,
                           show_cluster_labels = FALSE,
                           show_legend = TRUE,
                           legend_position = "right",
                           cluster_palette = NULL) {

  selection_mode <- match.arg(selection_mode)

  required_res_cols <- c(ontology_col, id_col, term_col, padj_col, size_col)
  missing_res_cols <- setdiff(required_res_cols, colnames(res_df))

  if (length(missing_res_cols) > 0) {
    stop(
      "Missing columns in res_df: ",
      paste(missing_res_cols, collapse = ", ")
    )
  }

  required_members_cols <- c("ID", "pathway_genes")
  missing_members_cols <- setdiff(required_members_cols, colnames(members_df))

  if (length(missing_members_cols) > 0) {
    stop(
      "Missing columns in members_df: ",
      paste(missing_members_cols, collapse = ", ")
    )
  }

  selected <- select_cerno_terms_for_tree(
    df = res_df,
    ontology = ontology,
    ontology_col = ontology_col,
    padj_col = padj_col,
    top_n = top_n,
    padj_cutoff = padj_cutoff,
    selection_mode = selection_mode
  )

  if (nrow(selected) < 3) {
    stop(
      "Too few terms for the tree plot in ontology ",
      ontology,
      ". At least three terms are required."
    )
  }

  dat <- selected |>
    dplyr::left_join(
      members_df,
      by = stats::setNames("ID", id_col)
    ) |>
    dplyr::filter(
      !purrr::map_lgl(pathway_genes, ~ is.null(.x) || length(.x) == 0)
    ) |>
    dplyr::mutate(
      tip_label = clean_pathway_name(
        .data[[term_col]],
        split_after_words = split_after_words
      ),
      phy_label = make.unique(
        stringr::str_replace_all(tip_label, "\n", " ")
      ),
      neglog10_padj = -log10(
        pmax(.data[[padj_col]], .Machine$double.xmin)
      )
    )

  if (nrow(dat) < 3) {
    stop("Fewer than three terms remain after joining members_df.")
  }

  n <- nrow(dat)

  sim_mat <- matrix(0, nrow = n, ncol = n)

  for (i in seq_len(n)) {
    for (j in i:n) {
      sim_ij <- jaccard_index(
        dat$pathway_genes[[i]],
        dat$pathway_genes[[j]]
      )
      sim_mat[i, j] <- sim_ij
      sim_mat[j, i] <- sim_ij
    }
  }

  rownames(sim_mat) <- dat$phy_label
  colnames(sim_mat) <- dat$phy_label

  dist_mat <- as.dist(1 - sim_mat)
  hc <- hclust(dist_mat, method = cluster_method)

  nCluster <- min(nCluster, nrow(dat) - 1)
  nCluster <- max(nCluster, 2)

  cluster_raw <- cutree(hc, k = nCluster)
  cluster_id <- paste0("c", cluster_raw)
  names(cluster_id) <- names(cluster_raw)

  phy <- ape::as.phylo(hc)

  tip_data <- dat |>
    dplyr::mutate(
      cluster = unname(cluster_id[phy_label])
    ) |>
    dplyr::transmute(
      label = phy_label,
      tip_label = tip_label,
      cluster = cluster,
      adj_pval = .data[[padj_col]],
      neglog10_padj = neglog10_padj,
      size_value = .data[[size_col]]
    )

  cluster_labels <- tip_data |>
    dplyr::group_by(cluster) |>
    dplyr::summarise(
      cluster_name = make_cluster_label(
        tip_label,
        max_words = cluster_label_words
      ),
      n_terms = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      cluster_name = make.unique(cluster_name)
    )

  cluster_label_map <- stats::setNames(
    cluster_labels$cluster_name,
    cluster_labels$cluster
  )

  tip_data <- tip_data |>
    dplyr::mutate(
      cluster_name = unname(cluster_label_map[cluster])
    )

  tip_groups <- split(tip_data$label, tip_data$cluster)

  clade_nodes <- purrr::map_dfr(
    names(tip_groups),
    function(cl) {
      tips <- tip_groups[[cl]]
      tip_idx <- match(tips, phy$tip.label)

      node <- if (length(tip_idx) == 1) {
        tip_idx
      } else {
        ape::getMRCA(phy, tip_idx)
      }

      tibble::tibble(
        cluster = cl,
        node = node,
        cluster_name = unname(cluster_label_map[cl])
      )
    }
  )

  p0 <- ggtree::ggtree(
    phy,
    layout = "rectangular",
    branch.length = branch_length,
    linewidth = branch_size
  )

  p0$data$x <- p0$data$x * tree_scale

  tree_max_x <- max(p0$data$x, na.rm = TRUE)

  tip_y <- p0$data |>
    dplyr::filter(isTip) |>
    dplyr::select(label, y)

  cluster_label_pos <- tip_data |>
    dplyr::left_join(tip_y, by = "label") |>
    dplyr::group_by(cluster, cluster_name) |>
    dplyr::summarise(
      y = mean(range(y, na.rm = TRUE)),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      x = tree_max_x + group_label_offset
    )

  x_right <- tree_max_x + tip_label_offset + tip_label_space + right_padding

  if (show_cluster_labels) {
    x_right <- max(
      x_right,
      max(cluster_label_pos$x, na.rm = TRUE) + right_padding
    )
  }

  cluster_levels <- sort(unique(tip_data$cluster_name))

  if (is.null(cluster_palette)) {
    cluster_palette <- default_cluster_palette(length(cluster_levels))
  }

  cluster_palette <- stats::setNames(
    cluster_palette[seq_along(cluster_levels)],
    cluster_levels
  )

  p <- p0 %<+% tip_data +
    ggtree::geom_hilight(
      data = clade_nodes,
      aes(
        node = node,
        fill = cluster_name
      ),
      alpha = clade_alpha,
      extend = hilight_extend,
      show.legend = show_legend
    ) +
    ggtree::geom_tippoint(
      aes(
        size = size_value,
        color = cluster_name
      ),
      alpha = 0.98,
      show.legend = TRUE
    ) +
    ggtree::geom_tiplab(
      aes(
        label = tip_label,
        color = cluster_name
      ),
      size = tip_label_size,
      offset = tip_label_offset,
      align = TRUE,
      linesize = 0.10,
      show.legend = FALSE
    ) +
    scale_fill_manual(
      values = cluster_palette,
      name = "Functional cluster"
    ) +
    scale_color_manual(
      values = cluster_palette,
      guide = "none"
    ) +
    scale_size_continuous(
      range = point_size_range,
      name = size_col
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(alpha = 0.65),
        order = 1
      ),
      size = guide_legend(order = 2)
    ) +
    coord_cartesian(
      xlim = c(0, x_right),
      clip = "off"
    ) +
    labs(
      title = paste0("CERNO GO:", ontology, " treeplot")
    ) +
    ggtree::theme_tree() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = if (show_legend) legend_position else "none",
      legend.justification = "center",
      legend.box = "vertical",
      legend.title = element_text(face = "bold", size = 9),
      legend.text = element_text(size = 8),
      legend.key.size = grid::unit(0.45, "cm"),
      legend.margin = margin(0, 0, 0, 8),
      legend.box.margin = margin(0, 0, 0, 14),
      plot.margin = margin(8, 8, 8, 8)
    )

  if (show_cluster_labels) {
    p <- p +
      geom_label(
        data = cluster_label_pos,
        aes(
          x = x,
          y = y,
          label = cluster_name,
          color = cluster_name
        ),
        inherit.aes = FALSE,
        hjust = 0,
        fontface = "bold",
        size = clade_label_size,
        label.size = 0,
        fill = "white",
        alpha = 0.92,
        show.legend = FALSE
      )
  }

  p
}
