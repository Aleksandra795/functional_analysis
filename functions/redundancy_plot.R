# GO semantic redundancy analysis and network plot

cerno_simplify_go_terms <- function(res_df,
                                    go_term_map,
                                    ontology = c("BP", "CC", "MF"),
                                    ontology_col = "ONTOLOGY",
                                    id_col = "ID",
                                    term_col = "Title",
                                    padj_col = "adj.P.Val",
                                    auc_col = "AUC",
                                    n_genes_col = "n_pathway_total",
                                    padj_cutoff = 0.05,
                                    top_n = NULL,
                                    selection_mode = c(
                                      "significant_top_n",
                                      "significant_all",
                                      "top_n",
                                      "significant_or_top_n"
                                    ),
                                    similarity_cutoff = 0.7,
                                    semantic_measure = "Wang",
                                    split_after_words = 7) {
  
  ontology <- match.arg(ontology)
  selection_mode <- match.arg(selection_mode)
  
  required_cols <- c(ontology_col, id_col, term_col, padj_col, auc_col)
  missing_cols <- setdiff(required_cols, colnames(res_df))
  
  if (length(missing_cols) > 0) {
    stop(
      "Missing columns in res_df: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  go_map <- go_term_map |>
    dplyr::transmute(
      map_id = ID,
      GOID = GOID,
      ONTOLOGY_map = ONTOLOGY_msigdbr
    ) |>
    dplyr::distinct(map_id, .keep_all = TRUE)
  
  x <- res_df |>
    dplyr::filter(.data[[ontology_col]] == ontology) |>
    dplyr::left_join(
      go_map,
      by = stats::setNames("map_id", id_col)
    ) |>
    dplyr::mutate(
      GOID = dplyr::case_when(
        grepl("^GO:", as.character(.data[[id_col]])) ~ as.character(.data[[id_col]]),
        !is.na(GOID) ~ GOID,
        TRUE ~ NA_character_
      )
    ) |>
    dplyr::filter(
      !is.na(GOID),
      !is.na(.data[[padj_col]]),
      !is.na(.data[[auc_col]])
    ) |>
    dplyr::arrange(.data[[padj_col]], dplyr::desc(.data[[auc_col]])) |>
    dplyr::distinct(GOID, .keep_all = TRUE) |>
    dplyr::mutate(
      pathway_clean = clean_pathway_name(
        .data[[term_col]],
        split_after_words = split_after_words
      ),
      neglog10_padj = -log10(
        pmax(.data[[padj_col]], .Machine$double.xmin)
      ),
      AUC_value = .data[[auc_col]],
      n_genes_value = if (n_genes_col %in% colnames(res_df)) {
        .data[[n_genes_col]]
      } else {
        NA_real_
      }
    )
  
  if (selection_mode == "significant_top_n") {
    
    if (is.null(top_n)) {
      stop("Set top_n when selection_mode = 'significant_top_n'.")
    }
    
    x <- x |>
      dplyr::filter(.data[[padj_col]] < padj_cutoff) |>
      dplyr::slice_head(n = top_n)
    
  } else if (selection_mode == "significant_all") {
    
    x <- x |>
      dplyr::filter(.data[[padj_col]] < padj_cutoff)
    
  } else if (selection_mode == "top_n") {
    
    if (is.null(top_n)) {
      stop("Set top_n when selection_mode = 'top_n'.")
    }
    
    x <- x |>
      dplyr::slice_head(n = top_n)
    
  } else if (selection_mode == "significant_or_top_n") {
    
    if (is.null(top_n)) {
      stop("Set top_n when selection_mode = 'significant_or_top_n'.")
    }
    
    sig <- x |>
      dplyr::filter(.data[[padj_col]] < padj_cutoff)
    
    if (nrow(sig) >= 2) {
      x <- sig |>
        dplyr::slice_head(n = top_n)
    } else {
      x <- x |>
        dplyr::slice_head(n = top_n)
    }
  }
  
  if (nrow(x) == 0) {
    stop("No GO terms remain after filtering.")
  }
  
  if (nrow(x) == 1) {
    
    x_out <- x |>
      dplyr::mutate(
        redundancy_cluster = "R01",
        cluster_size = 1L,
        is_representative = TRUE,
        representative_GOID = GOID,
        representative_label = pathway_clean,
        redundant_to = NA_character_
      )
    
    return(
      list(
        all_terms = x_out,
        representatives = x_out,
        removed = x_out[0, ],
        similarity_matrix = matrix(1, nrow = 1, ncol = 1, dimnames = list(x$GOID, x$GOID)),
        similarity_cutoff = similarity_cutoff,
        ontology = ontology,
        semantic_measure = semantic_measure,
        summary = tibble::tibble(
          ontology = ontology,
          n_input_terms = 1L,
          n_representatives = 1L,
          n_removed = 0L,
          similarity_cutoff = similarity_cutoff,
          semantic_measure = semantic_measure
        )
      )
    )
  }
  
  sim_mat <- compute_go_semantic_similarity(
    go_ids = x$GOID,
    ontology = ontology,
    measure = semantic_measure
  )
  
  adjacency <- sim_mat >= similarity_cutoff
  diag(adjacency) <- FALSE
  
  g <- igraph::graph_from_adjacency_matrix(
    adjacency,
    mode = "undirected",
    diag = FALSE
  )
  
  comp <- igraph::components(g)$membership
  
  cluster_id <- sprintf(
    "R%02d",
    as.integer(factor(comp[x$GOID]))
  )
  
  x_annotated <- x |>
    dplyr::mutate(
      redundancy_cluster = cluster_id
    ) |>
    dplyr::group_by(redundancy_cluster) |>
    dplyr::arrange(
      .data[[padj_col]],
      dplyr::desc(.data[[auc_col]]),
      .by_group = TRUE
    ) |>
    dplyr::mutate(
      cluster_size = dplyr::n(),
      is_representative = dplyr::row_number() == 1,
      representative_GOID = dplyr::first(GOID),
      representative_label = dplyr::first(pathway_clean),
      redundant_to = dplyr::if_else(
        is_representative,
        NA_character_,
        representative_label
      )
    ) |>
    dplyr::ungroup()
  
  representatives <- x_annotated |>
    dplyr::filter(is_representative) |>
    dplyr::arrange(.data[[padj_col]])
  
  removed <- x_annotated |>
    dplyr::filter(!is_representative) |>
    dplyr::arrange(redundant_to, .data[[padj_col]])
  
  list(
    all_terms = x_annotated,
    representatives = representatives,
    removed = removed,
    similarity_matrix = sim_mat,
    similarity_cutoff = similarity_cutoff,
    ontology = ontology,
    semantic_measure = semantic_measure,
    summary = tibble::tibble(
      ontology = ontology,
      n_input_terms = nrow(x_annotated),
      n_representatives = nrow(representatives),
      n_removed = nrow(removed),
      similarity_cutoff = similarity_cutoff,
      semantic_measure = semantic_measure
    )
  )
}

compute_go_semantic_similarity <- function(go_ids,
                                           ontology = c("BP", "CC", "MF"),
                                           measure = "Wang") {
  
  ontology <- match.arg(ontology)
  go_ids <- unique(as.character(go_ids))
  
  compute_ic <- !(measure %in% c("Wang"))
  
  sem_data <- GOSemSim::godata(
    annoDb = "org.Hs.eg.db",
    ont = ontology,
    computeIC = compute_ic
  )
  
  n <- length(go_ids)
  
  sim_mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(go_ids, go_ids)
  )
  
  diag(sim_mat) <- 1
  
  if (n < 2) {
    return(sim_mat)
  }
  
  for (i in seq_len(n)) {
    for (j in i:n) {
      
      if (i == j) {
        sim_ij <- 1
      } else {
        sim_ij <- GOSemSim::goSim(
          go_ids[i],
          go_ids[j],
          semData = sem_data,
          measure = measure
        )
        
        sim_ij <- as.numeric(sim_ij)[1]
        
        if (is.na(sim_ij) || is.nan(sim_ij)) {
          sim_ij <- 0
        }
        
        sim_ij <- max(0, min(1, sim_ij))
      }
      
      sim_mat[i, j] <- sim_ij
      sim_mat[j, i] <- sim_ij
    }
  }
  
  sim_mat
}

cerno_go_redundancy_emap <- function(simplify_result,
                                     label = c("representatives", "all", "none"),
                                     similarity_cutoff = NULL,
                                     layout = "fr",
                                     label_size = 3,
                                     show_removed_labels = TRUE,
                                     removed_label_color = "grey50",
                                     removed_label_alpha = 0.55,
                                     removed_label_size_factor = 0.75,
                                     edge_alpha = 0.45,
                                     edge_width_range = c(0.3, 2.8),
                                     node_size_range = c(3, 9),
                                     box_padding=0.7,
                                     title = NULL) {
  
  label <- match.arg(label)
  
  terms <- simplify_result$all_terms
  sim_mat <- simplify_result$similarity_matrix
  
  if (is.null(similarity_cutoff)) {
    similarity_cutoff <- simplify_result$similarity_cutoff
  }
  
  idx <- which(
    sim_mat >= similarity_cutoff & upper.tri(sim_mat),
    arr.ind = TRUE
  )
  
  if (nrow(idx) > 0) {
    edges <- tibble::tibble(
      from = rownames(sim_mat)[idx[, 1]],
      to = colnames(sim_mat)[idx[, 2]],
      similarity = as.numeric(sim_mat[idx])
    )
  } else {
    edges <- tibble::tibble(
      from = character(),
      to = character(),
      similarity = numeric()
    )
  }
  
  vertices <- terms |>
    dplyr::transmute(
      name = GOID,
      pathway_clean = pathway_clean,
      redundancy_cluster = redundancy_cluster,
      cluster_size = cluster_size,
      is_representative = is_representative,
      neglog10_padj = neglog10_padj,
      AUC_value = AUC_value
    ) |>
    dplyr::mutate(
      annotation_type = dplyr::case_when(
        label %in% c("representatives", "all") & is_representative ~ "representative",
        show_removed_labels & label != "none" & !is_representative ~ "removed",
        TRUE ~ "none"
      ),
      annotation_label = dplyr::case_when(
        annotation_type != "none" ~ pathway_clean,
        TRUE ~ ""
      ),
      annotation_fontface = dplyr::case_when(
        annotation_type == "representative" ~ "bold",
        TRUE ~ "plain"
      )
    )
  
  annotation_colors <- c(
    representative = "black",
    removed = grDevices::adjustcolor(
      removed_label_color,
      alpha.f = removed_label_alpha
    )
  )
  
  if (nrow(edges) > 0) {
    
    similarity_breaks <- sort(unique(c(
      similarity_cutoff,
      0.80,
      0.90,
      Inf
    )))
    
    if (length(similarity_breaks) < 2) {
      similarity_breaks <- c(similarity_cutoff, Inf)
    }
    
    fmt <- function(x) {
      sprintf("%.2f", x)
    }
    
    similarity_labels <- vapply(
      seq_len(length(similarity_breaks) - 1),
      function(i) {
        left <- similarity_breaks[i]
        right <- similarity_breaks[i + 1]
        
        if (is.infinite(right)) {
          paste0("\u2265", fmt(left))
        } else {
          paste0(fmt(left), "\u2013<", fmt(right))
        }
      },
      character(1)
    )
    
    edges <- edges |>
      dplyr::mutate(
        similarity_class = cut(
          similarity,
          breaks = similarity_breaks,
          labels = similarity_labels,
          include.lowest = TRUE,
          right = FALSE
        ),
        similarity_class = factor(
          similarity_class,
          levels = similarity_labels
        )
      )
    
    present_similarity_levels <- levels(droplevels(edges$similarity_class))
    
    edge_width_values <- stats::setNames(
      seq(
        edge_width_range[1],
        edge_width_range[2],
        length.out = length(present_similarity_levels)
      ),
      present_similarity_levels
    )
    
    edge_legend_df <- tibble::tibble(
      x = 0,
      xend = 0,
      y = 0,
      yend = 0,
      similarity_class = factor(
        present_similarity_levels,
        levels = present_similarity_levels
      )
    )
    
  } else {
    
    present_similarity_levels <- character(0)
    edge_width_values <- numeric(0)
    
    edge_legend_df <- tibble::tibble(
      x = numeric(),
      xend = numeric(),
      y = numeric(),
      yend = numeric(),
      similarity_class = factor()
    )
  }
  
  graph <- igraph::graph_from_data_frame(
    d = edges |>
      dplyr::select(from, to, similarity, dplyr::any_of("similarity_class")),
    directed = FALSE,
    vertices = vertices
  )
  
  if (is.null(title)) {
    title <- paste0(
      "CERNO GO:",
      simplify_result$ontology,
      " redundancy network"
    )
  }
  
  p <- ggraph::ggraph(graph, layout = layout)
  
  if (nrow(edges) > 0) {
    p <- p +
      ggraph::geom_edge_link(
        ggplot2::aes(edge_width = similarity_class),
        color = "grey55",
        alpha = edge_alpha,
        show.legend = FALSE
      ) +
      ggraph::scale_edge_width_manual(
        values = edge_width_values,
        guide = "none"
      ) +
      ggplot2::geom_segment(
        data = edge_legend_df,
        ggplot2::aes(
          x = x,
          xend = xend,
          y = y,
          yend = yend,
          linewidth = similarity_class
        ),
        inherit.aes = FALSE,
        color = "grey55",
        alpha = 0,
        show.legend = TRUE,
        key_glyph = ggplot2::draw_key_path
      ) +
      ggplot2::scale_linewidth_manual(
        values = edge_width_values,
        name = "Semantic similarity"
      )
  }
  
  p +
    ggraph::geom_node_point(
      ggplot2::aes(
        size = AUC_value,
        fill = neglog10_padj,
        alpha = is_representative
      ),
      shape = 21,
      color = "black",
      stroke = 0.35,
      show.legend = c(
        size = TRUE,
        fill = TRUE,
        alpha = FALSE
      ),
      key_glyph = ggplot2::draw_key_point
    ) +
    ggraph::geom_node_text(
      data = function(x) dplyr::filter(x, annotation_label != ""),
      ggplot2::aes(
        label = annotation_label,
        color = annotation_type,
        fontface = annotation_fontface
      ),
      repel = TRUE,
      size = label_size,
      box.padding = grid::unit(box_padding, "lines"),
      point.padding = grid::unit(0.15, "lines"),
      force = 3,
      force_pull = 0.8,
      max.time = 15,
      max.iter = 120000,
      min.segment.length = grid::unit(1.5, "lines"),
      segment.color = "grey70",
      segment.alpha = 0.35,
      segment.size = 0.15,
      show.legend = FALSE,
      max.overlaps = Inf
    ) +
    ggplot2::scale_color_manual(
      values = annotation_colors,
      guide = "none"
    ) +
    ggplot2::scale_alpha_manual(
      values = c(
        `TRUE` = 1,
        `FALSE` = 0.30
      ),
      guide = "none"
    ) +
    ggplot2::scale_size_continuous(
      range = node_size_range,
      name = "AUC"
    ) +
    ggplot2::scale_fill_gradient(
      low = "#2166AC",
      high = "#B2182B",
      name = expression(-log[10]("FDR"))
    ) +
    ggplot2::guides(
      linewidth = ggplot2::guide_legend(
        title = "Semantic similarity",
        order = 1,
        override.aes = list(
          alpha = 1,
          colour = "grey55"
        )
      ),
      fill = ggplot2::guide_colorbar(
        title = expression(-log[10]("FDR")),
        order = 2
      ),
      size = ggplot2::guide_legend(
        title = "AUC",
        order = 3,
        override.aes = list(
          shape = 21,
          colour = "black",
          fill = "grey80",
          alpha = 1,
          stroke = 0.35
        )
      )
    ) +
    ggplot2::labs(
      title = title
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "bold"),
      legend.text = ggplot2::element_text(size = 8),
      plot.margin = ggplot2::margin(10, 20, 10, 20)
    )
}
