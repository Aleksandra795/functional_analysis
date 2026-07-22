# Shared helpers

to_num <- function(x) {
  if (is.numeric(x)) {
    return(x)
  }
  as.numeric(gsub(",", ".", as.character(x)))
}

read_go_cerno_excel <- function(file,
                                ontology,
                                sheet = 1) {
  
  x <- readxl::read_excel(
    path = file,
    sheet = sheet
  )
  
  numeric_cols <- intersect(
    c(
      "AUC",
      "P.Value",
      "adj.P.Val",
      "n_overlap",
      "n_pathway_total",
      "coverage_pathway"
    ),
    names(x)
  )
  
  x <- x %>%
    mutate(
      ONTOLOGY = ontology,
      source_file = basename(file),
      .before = 1
    ) %>%
    mutate(
      across(all_of(numeric_cols), to_num)
    )
  x
}

load_msigdb_collections <- function(species = "Homo sapiens") {
  list(
    BP = msigdbr::msigdbr(
      species = species,
      collection = "C5",
      subcollection = "GO:BP"
    ),
    CC = msigdbr::msigdbr(
      species = species,
      collection = "C5",
      subcollection = "GO:CC"
    ),
    MF = msigdbr::msigdbr(
      species = species,
      collection = "C5",
      subcollection = "GO:MF"
    ),
    REACTOME = msigdbr::msigdbr(
      species = species,
      collection = "C2",
      subcollection = "CP:REACTOME"
    )
  )
}

normalize_pathway_key <- function(x) {
  x <- as.character(x)
  x <- toupper(x)
  x <- gsub("[^A-Z0-9]+", "_", x)
  gsub("^_|_$", "", x)
}

canonicalize_result_ids <- function(results,
                                    msig_df,
                                    id_col = "ID",
                                    term_col = "Title") {
  alias_map <- dplyr::bind_rows(
    dplyr::transmute(msig_df, alias = gs_id, canonical_ID = gs_id),
    dplyr::transmute(msig_df, alias = gs_name, canonical_ID = gs_id),
    dplyr::transmute(msig_df, alias = gs_exact_source, canonical_ID = gs_id)
  ) |>
    dplyr::filter(!is.na(alias), alias != "") |>
    dplyr::mutate(alias_key = normalize_pathway_key(alias)) |>
    dplyr::distinct(alias_key, .keep_all = TRUE)
  
  x <- results |>
    dplyr::mutate(
      input_ID = as.character(.data[[id_col]]),
      id_key = normalize_pathway_key(.data[[id_col]]),
      term_key = normalize_pathway_key(.data[[term_col]])
    ) |>
    dplyr::left_join(
      dplyr::select(alias_map, id_key = alias_key, canonical_from_id = canonical_ID),
      by = "id_key"
    ) |>
    dplyr::left_join(
      dplyr::select(alias_map, term_key = alias_key, canonical_from_term = canonical_ID),
      by = "term_key"
    ) |>
    dplyr::mutate(
      "{id_col}" := dplyr::coalesce(
        canonical_from_id,
        canonical_from_term,
        as.character(.data[[id_col]])
      )
    ) |>
    dplyr::select(-id_key, -term_key, -canonical_from_id, -canonical_from_term)
  
  x
}

build_pathway_members <- function(msig_df) {
  msig_df |>
    dplyr::filter(!is.na(gs_id), !is.na(gene_symbol)) |>
    dplyr::distinct(gs_id, gs_name, gs_exact_source, gene_symbol) |>
    dplyr::group_by(gs_id, gs_name, gs_exact_source) |>
    dplyr::summarise(
      pathway_genes = list(sort(unique(gene_symbol))),
      n_pathway_total_msigdbr = dplyr::n_distinct(gene_symbol),
      .groups = "drop"
    ) |>
    dplyr::transmute(
      ID = gs_id,
      Title_msigdbr = gs_name,
      source_id = gs_exact_source,
      pathway_genes,
      n_pathway_total_msigdbr
    )
}

build_go_term_map <- function(msig_go_df) {
  msig_go_df |>
    dplyr::distinct(gs_id, gs_name, gs_exact_source, gs_subcollection) |>
    dplyr::transmute(
      ID = gs_id,
      GOID = gs_exact_source,
      Title_msigdbr = gs_name,
      ONTOLOGY_msigdbr = dplyr::case_when(
        gs_subcollection == "GO:BP" ~ "BP",
        gs_subcollection == "GO:CC" ~ "CC",
        gs_subcollection == "GO:MF" ~ "MF",
        TRUE ~ gs_subcollection
      )
    ) |>
    dplyr::filter(grepl("^GO:", GOID))
}


filter_pathways_and_recalculate_fdr <- function(
    df,
    min_pathway_size = 10,
    max_pathway_size = 1099,
    size_col = "n_pathway_total",
    pval_col = "P.Value",
    padj_col = "adj.P.Val",
    original_padj_col = "adj.P.Val_original",
    method = "fdr"
) {
  x <- df
  
  if (!original_padj_col %in% colnames(x)) {
    x[[original_padj_col]] <- x[[padj_col]]
  }
  
  x <- x |>
    dplyr::filter(
      .data[[size_col]] >= min_pathway_size,
      .data[[size_col]] <= max_pathway_size
    )
  
  x[[padj_col]] <- stats::p.adjust(
    x[[pval_col]],
    method = method
  )
  
  x |>
    dplyr::arrange(
      .data[[padj_col]],
      .data[[pval_col]]
    )
}


clean_pathway_name <- function(x, split_after_words = 15) {

  x <- as.character(x)

  x <- gsub(
    "^(GOBP_|GOCC_|GOMF_|REACTOME_|KEGG_|HALLMARK_|BIOCARTA_|PID_|WP_)",
    "",
    x
  )

  x <- gsub("_", " ", x)

  x <- tolower(x)

  x <- stringr::str_squish(x)

  x <- vapply(
    x,
    function(term) {
      words <- unlist(strsplit(term, "\\s+"))
      n_words <- length(words)

      if (n_words > split_after_words) {
        cut_point <- ceiling(n_words / 2)
        paste0(
          paste(words[1:cut_point], collapse = " "),
          "\n",
          paste(words[(cut_point + 1):n_words], collapse = " ")
        )
      } else {
        term
      }
    },
    character(1)
  )

  x
}

highlight_pathway_labels <- function(labels,
                                     highlight1 = NULL,
                                     highlight2 = NULL,
                                     color1 = "#009E73",
                                     color2 = "#0072B2") {

  out <- labels

  if (!is.null(highlight1)) {
    highlight1_clean <- clean_pathway_name(highlight1)
    idx1 <- labels %in% highlight1_clean

    out[idx1] <- paste0(
      "<span style='color:", color1, "; font-weight:700;'>",
      labels[idx1],
      "</span>"
    )
  }

  if (!is.null(highlight2)) {
    highlight2_clean <- clean_pathway_name(highlight2)
    idx2 <- labels %in% highlight2_clean

    out[idx2] <- paste0(
      "<span style='color:", color2, "; font-weight:700;'>",
      labels[idx2],
      "</span>"
    )
  }

  out
}

jaccard_index <- function(a, b) {
  u <- union(a, b)
  if (length(u) == 0) return(0)
  length(intersect(a, b)) / length(u)
}

build_jaccard_matrix <- function(gene_sets, labels) {

  n <- length(gene_sets)
  sim_mat <- matrix(0, nrow = n, ncol = n)

  for (i in seq_len(n)) {
    for (j in i:n) {
      s <- jaccard_index(gene_sets[[i]], gene_sets[[j]])
      sim_mat[i, j] <- s
      sim_mat[j, i] <- s
    }
  }

  rownames(sim_mat) <- labels
  colnames(sim_mat) <- labels

  sim_mat
}

make_cluster_label <- function(terms,
                               max_words = 3,
                               stop_words = NULL) {

  if (is.null(stop_words)) {
    stop_words <- c(
      "go", "biological", "process", "molecular", "function",
      "cellular", "component", "regulation", "positive", "negative",
      "response", "activity", "binding", "protein",
      "cell", "cells", "of", "to", "in", "by", "via", "and",
      "or", "from", "with", "dependent", "mediated", "pathway"
    )
  }

  words <- terms |>
    str_replace_all("\n", " ") |>
    str_to_lower() |>
    str_split("\\s+") |>
    unlist()

  words <- gsub("[^a-z0-9]", "", words)
  words <- words[nchar(words) > 2]
  words <- words[!(words %in% stop_words)]

  if (length(words) == 0) return("functional cluster")

  tab <- sort(table(words), decreasing = TRUE)
  top_words <- names(tab)[seq_len(min(max_words, length(tab)))]

  paste(top_words, collapse = " / ")
}

default_cluster_palette <- function(n) {
  base_cols <- c(
    "#1B9E77", "#D95F02", "#7570B3", "#E7298A",
    "#66A61E", "#E6AB02", "#A6761D", "#1F78B4",
    "#B2DF8A", "#FB9A99", "#6A3D9A", "#B15928"
  )

  rep(base_cols, length.out = n)
}
