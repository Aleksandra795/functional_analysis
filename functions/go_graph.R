# GO hierarchy graph

get_go_parent_map <- function(ontology = c("BP", "CC", "MF")) {

  ontology <- match.arg(ontology)

  parent_env <- switch(
    ontology,
    BP = GO.db::GOBPPARENTS,
    CC = GO.db::GOCCPARENTS,
    MF = GO.db::GOMFPARENTS
  )

  AnnotationDbi::as.list(parent_env)
}

normalize_go_relation <- function(x) {

  x <- as.character(x)
  x <- stringr::str_to_lower(x)
  x <- gsub("[^a-z0-9]+", "_", x)
  x <- gsub("^_|_$", "", x)

  dplyr::case_when(
    x %in% c("isa", "is_a") ~ "is_a",
    x %in% c("partof", "part_of") ~ "part_of",
    x %in% c("regulates", "regulates_process") ~ "regulates",
    x %in% c("positively_regulates", "positively_regulates_process") ~ "positively_regulates",
    x %in% c("negatively_regulates", "negatively_regulates_process") ~ "negatively_regulates",
    TRUE ~ "parent"
  )
}

get_go_ancestor_edges <- function(go_ids,
                                  ontology = c("BP", "CC", "MF"),
                                  max_depth = Inf,
                                  include_roots = TRUE) {

  ontology <- match.arg(ontology)
  parent_map <- get_go_parent_map(ontology)

  root_id <- switch(
    ontology,
    BP = "GO:0008150",
    CC = "GO:0005575",
    MF = "GO:0003674"
  )

  visited <- character(0)
  frontier <- unique(go_ids)

  edges <- tibble::tibble(
    from = character(),
    to = character(),
    relation = character(),
    depth = integer()
  )

  depth <- 0

  while (length(frontier) > 0 && depth < max_depth) {

    depth <- depth + 1
    next_frontier <- character(0)

    for (child in frontier) {

      if (child %in% visited) next
      visited <- c(visited, child)

      parents <- parent_map[[child]]
      if (is.null(parents) || length(parents) == 0) next

      parent_values <- as.character(parents)
      parent_names <- names(parents)

      if (any(grepl("^GO:", parent_values))) {
        parent_ids <- parent_values

        if (!is.null(parent_names) && length(parent_names) == length(parent_ids)) {
          relation_raw <- parent_names
        } else {
          relation_raw <- rep("parent", length(parent_ids))
        }

      } else if (!is.null(parent_names) && any(grepl("^GO:", parent_names))) {
        parent_ids <- parent_names
        relation_raw <- parent_values

      } else {
        next
      }

      relation_types <- normalize_go_relation(relation_raw)

      keep <- grepl("^GO:", parent_ids)

      if (!include_roots) {
        keep <- keep & parent_ids != root_id
      }

      parent_ids <- parent_ids[keep]
      relation_types <- relation_types[keep]

      if (length(parent_ids) == 0) next

      edges <- dplyr::bind_rows(
        edges,
        tibble::tibble(
          from = parent_ids,
          to = child,
          relation = relation_types,
          depth = depth
        )
      )

      next_frontier <- c(next_frontier, parent_ids)
    }

    frontier <- unique(next_frontier)
  }

  edges |>
    dplyr::distinct(from, to, relation, .keep_all = TRUE)
}

get_go_terms_metadata <- function(go_ids) {

  meta <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(go_ids),
    keytype = "GOID",
    columns = c("GOID", "TERM", "ONTOLOGY")
  )

  meta |>
    distinct(GOID, .keep_all = TRUE)
}

select_cerno_go_terms_for_graph <- function(res_df,
                                            go_term_map,
                                            ontology = c("BP", "CC", "MF"),
                                            ontology_col = "ONTOLOGY",
                                            id_col = "ID",
                                            term_col = "Title",
                                            padj_col = "adj.P.Val",
                                            auc_col = "AUC",
                                            n_genes_col = "n_pathway_total",
                                            top_n = 15,
                                            padj_cutoff = 0.05,
                                            selection_mode = c(
                                              "significant_top_n",
                                              "significant_all",
                                              "top_n",
                                              "significant_or_top_n"
                                            )) {

  ontology <- match.arg(ontology)
  selection_mode <- match.arg(selection_mode)

  required_cols <- c(ontology_col, id_col, term_col, padj_col, auc_col)
  missing_cols <- setdiff(required_cols, colnames(res_df))

  if (length(missing_cols) > 0) {
    stop("Missing columns in res_df: ", paste(missing_cols, collapse = ", "))
  }

  x <- res_df |>
    filter(.data[[ontology_col]] == ontology) |>
    left_join(go_term_map, by = setNames("ID", id_col)) |>
    mutate(
      GOID = dplyr::case_when(
        grepl("^GO:", as.character(.data[[id_col]])) ~ as.character(.data[[id_col]]),
        !is.na(GOID) ~ GOID,
        TRUE ~ NA_character_
      )
    ) |>
    filter(!is.na(GOID), !is.na(.data[[padj_col]])) |>
    arrange(.data[[padj_col]])

  if (selection_mode == "significant_top_n") {

    x <- x |>
      filter(.data[[padj_col]] < padj_cutoff) |>
      slice_head(n = top_n)

  } else if (selection_mode == "significant_all") {

    x <- x |>
      filter(.data[[padj_col]] < padj_cutoff)

  } else if (selection_mode == "top_n") {

    x <- x |>
      slice_head(n = top_n)

  } else if (selection_mode == "significant_or_top_n") {

    sig <- x |>
      filter(.data[[padj_col]] < padj_cutoff)

    if (nrow(sig) >= 2) {
      x <- sig |>
        slice_head(n = top_n)
    } else {
      x <- x |>
        slice_head(n = top_n)
    }
  }

  if (nrow(x) == 0) {
    stop("No GO terms remain after filtering.")
  }

  x |>
    mutate(
      result_label = clean_pathway_name(.data[[term_col]], split_after_words = 6),
      neglog10_padj = -log10(pmax(.data[[padj_col]], .Machine$double.xmin)),
      AUC_value = .data[[auc_col]],
      n_genes_value = if (n_genes_col %in% colnames(x)) .data[[n_genes_col]] else NA_real_
    )
}

cerno_go_graph <- function(res_df,
                           go_term_map,
                           ontology = c("BP", "CC", "MF"),
                           ontology_col = "ONTOLOGY",
                           id_col = "ID",
                           term_col = "Title",
                           padj_col = "adj.P.Val",
                           auc_col = "AUC",
                           n_genes_col = "n_pathway_total",
                           top_n = 12,
                           padj_cutoff = 0.05,
                           selection_mode = c(
                             "significant_top_n",
                             "significant_all",
                             "top_n",
                             "significant_or_top_n"
                           ),
                           max_depth = 4,
                           include_roots = TRUE,
                           label_ancestors = FALSE,
                           label_size = 3,
                           node_size_range = c(3, 9),
                           edge_alpha = 0.90,
                           edge_width = 0.75,
                           ancestor_node_size = 2.6,
                           ancestor_node_fill = "grey80",
                           ancestor_node_color = "grey35",
                           layout = "sugiyama",
                           title = NULL,
                           show_components_in_title = FALSE) {

  ontology <- match.arg(ontology)
  selection_mode <- match.arg(selection_mode)

  selected_terms <- select_cerno_go_terms_for_graph(
    res_df = res_df,
    go_term_map = go_term_map,
    ontology = ontology,
    ontology_col = ontology_col,
    id_col = id_col,
    term_col = term_col,
    padj_col = padj_col,
    auc_col = auc_col,
    n_genes_col = n_genes_col,
    top_n = top_n,
    padj_cutoff = padj_cutoff,
    selection_mode = selection_mode
  )

  edges <- get_go_ancestor_edges(
    go_ids = selected_terms$GOID,
    ontology = ontology,
    max_depth = max_depth,
    include_roots = include_roots
  )

  if (nrow(edges) == 0) {
    stop("No GO edges could be built for the selected terms.")
  }

  edge_linetypes_all <- c(
    "is_a" = "solid",
    "part_of" = "dashed",
    "regulates" = "dotdash",
    "positively_regulates" = "twodash",
    "negatively_regulates" = "dotted",
    "parent" = "solid"
  )

  edge_colours_all <- c(
    "is_a" = "grey20",
    "part_of" = "black",
    "regulates" = "#542788",
    "positively_regulates" = "#006D2C",
    "negatively_regulates" = "#A63603",
    "parent" = "grey35"
  )

  edges <- edges |>
    dplyr::mutate(
      relation = as.character(relation),
      relation = dplyr::if_else(
        relation %in% names(edge_linetypes_all),
        relation,
        "parent"
      )
    )

  present_relations <- sort(unique(edges$relation))

  edge_linetypes <- edge_linetypes_all[present_relations]
  edge_colours <- edge_colours_all[present_relations]

  node_ids <- unique(c(edges$from, edges$to, selected_terms$GOID))

  go_meta <- get_go_terms_metadata(node_ids)

  nodes <- tibble::tibble(GOID = node_ids) |>
    dplyr::left_join(go_meta, by = "GOID") |>
    dplyr::left_join(
      selected_terms |>
        dplyr::transmute(
          GOID,
          result_title = .data[[term_col]],
          result_label,
          adj_pval = .data[[padj_col]],
          neglog10_padj,
          AUC_value,
          n_genes_value
        ),
      by = "GOID"
    ) |>
    dplyr::mutate(
      is_result = !is.na(adj_pval),
      node_label = dplyr::case_when(
        is_result ~ result_label,
        label_ancestors ~ clean_pathway_name(TERM, split_after_words = 6),
        TRUE ~ ""
      )
    )

  graph <- igraph::graph_from_data_frame(
    d = edges |>
      dplyr::transmute(
        from = from,
        to = to,
        relation = relation
      ),
    directed = TRUE,
    vertices = nodes |>
      dplyr::rename(name = GOID)
  )

  n_components <- igraph::components(
    igraph::as.undirected(graph, mode = "collapse")
  )$no

  if (is.null(title)) {
    title <- paste0("CERNO GO:", ontology, " ontology graph")

    if (show_components_in_title) {
      title <- paste0(title, "; components = ", n_components)
    }
  }

  p <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        edge_colour = relation,
        edge_linetype = relation
      ),
      alpha = edge_alpha,
      linewidth = edge_width,
      arrow = grid::arrow(
        length = grid::unit(2.2, "mm"),
        type = "closed"
      ),
      end_cap = ggraph::circle(2.8, "mm"),
      show.legend = c(
        edge_colour = TRUE,
        edge_linetype = TRUE
      )
    ) +
    ggraph::scale_edge_colour_manual(
      values = edge_colours,
      name = "GO relationship"
    ) +
    ggraph::scale_edge_linetype_manual(
      values = edge_linetypes,
      name = "GO relationship"
    ) +
    ggraph::geom_node_point(
      data = function(x) dplyr::filter(x, !is_result),
      ggplot2::aes(
        x = x,
        y = y
      ),
      inherit.aes = FALSE,
      shape = 21,
      color = ancestor_node_color,
      fill = ancestor_node_fill,
      size = ancestor_node_size,
      stroke = 0.35,
      alpha = 0.90,
      show.legend = FALSE
    ) +
    ggraph::geom_node_point(
      data = function(x) dplyr::filter(x, is_result),
      ggplot2::aes(
        x = x,
        y = y,
        size = AUC_value,
        fill = neglog10_padj
      ),
      inherit.aes = FALSE,
      shape = 21,
      color = "black",
      stroke = 0.35,
      alpha = 0.95,
      show.legend = c(
        size = TRUE,
        fill = TRUE
      )
    ) +
    ggraph::geom_node_label(
      data = function(x) dplyr::filter(x, node_label != ""),
      ggplot2::aes(
        x = x,
        y = y,
        label = node_label
      ),
      inherit.aes = FALSE,
      repel = TRUE,
      size = label_size,
      label.size = 0.15,
      fill = "white",
      alpha = 0.90,
      max.overlaps = Inf,
      show.legend = FALSE
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
      edge_colour = ggplot2::guide_legend(
        title = "GO relationship",
        order = 1,
        override.aes = list(
          alpha = 1,
          linewidth = 0.9
        )
      ),
      edge_linetype = ggplot2::guide_legend(
        title = "GO relationship",
        order = 1,
        override.aes = list(
          alpha = 1,
          linewidth = 0.9
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

  p
}
