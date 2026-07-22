# CERNO functional enrichment visualizations
# Edit only the settings in the section below before running the script.

# Required packages -----------------------------------------------------------
library(AnnotationDbi)
library(ape)
library(dplyr)
library(GO.db)
library(GOSemSim)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(ggtext)
library(ggtree)
library(igraph)
library(msigdbr)
library(org.Hs.eg.db)
library(purrr)
library(readr)
library(stringr)
library(tibble)

# Local functions -------------------------------------------------------------
function_files <- c(
  "utils.R",
  "dotplot.R",
  "manhattan_plot.R",
  "treeplot.R",
  "go_graph.R",
  "redundancy_plot.R"
)

invisible(lapply(
  file.path("functions", function_files),
  source
))

# User settings ---------------------------------------------------------------
input_dir <- "path/to/CERNO/results"
output_dir <- file.path(input_dir, "plots")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# File names can be changed independently of the directory.
go_files <- c(
  BP = "CERNO_GO_BP.xlsx",
  CC = "CERNO_GO_CC.xlsx",
  MF = "CERNO_GO_MF.xlsx"
)

# Set to NULL to skip Reactome visualizations.
reactome_file <- "CERNO_REACTOME.xlsx"

plot_ontology <- "MF"       # one of: BP, CC, MF
padj_cutoff <- 0.05
top_n <- 30

# Optional pathway highlighting in the dot plot.
highlight_group_1 <- NULL
highlight_group_2 <- NULL

# MSigDB pathway definitions --------------------------------------------------
collections <- load_msigdb_collections()
go_msig <- bind_rows(collections$BP, collections$CC, collections$MF)
go_members <- build_pathway_members(go_msig)
go_term_map <- build_go_term_map(go_msig)
reactome_members <- build_pathway_members(collections$REACTOME)

# Read and normalize CERNO results -------------------------------------------
go_results <- purrr::imap_dfr(
  go_files,
  function(file_name, ontology) {
    read_go_cerno_excel(
      file = file.path(input_dir, file_name),
      ontology = ontology
    )
  }
)

go_results <- purrr::map_dfr(
  valid_ontologies,
  function(ontology) {
    canonicalize_result_ids(
      results = dplyr::filter(go_results, ONTOLOGY == ontology),
      msig_df = collections[[ontology]]
    )
  }
)

reactome_results <- NULL
if (!is.null(reactome_file)) {
  reactome_results <- read_go_cerno_excel(
    file = file.path(input_dir, reactome_file),
    ontology = "REACTOME"
  )

  reactome_results <- canonicalize_result_ids(
    results = reactome_results,
    msig_df = collections$REACTOME
  )
}

# Dot plot --------------------------------------------------------------------
p_dot <- cerno_dotplot(
  df = dplyr::filter(go_results, ONTOLOGY == plot_ontology),
  top_n = top_n,
  padj_cutoff = padj_cutoff,
  highlight1 = highlight_group_1,
  highlight2 = highlight_group_2,
  point_size_range = c(1.5, 7)
)

ggsave(
  file.path(output_dir, paste0("cerno_dotplot_", tolower(plot_ontology), ".png")),
  p_dot,
  width = 10,
  height = 8,
  dpi = 300
)

# Manhattan-style plot across GO ontologies ----------------------------------
p_manhattan <- cerno_manhattan_go(
  df = go_results,
  padj_cutoff = padj_cutoff,
  label_top_n = 15,
  point_size_range = c(0.8, 4),
  seed = 123
)

ggsave(
  file.path(output_dir, "cerno_manhattan_go.png"),
  p_manhattan,
  width = 11,
  height = 7,
  dpi = 300
)

# Pathway-overlap tree ---------------------------------------------------------
p_tree <- cerno_treeplot(
  res_df = go_results,
  members_df = go_members,
  ontology = plot_ontology,
  top_n = top_n,
  padj_cutoff = padj_cutoff,
  selection_mode = "significant_or_top_n",
  nCluster = 5,
  show_cluster_labels = FALSE,
  show_legend = TRUE
)

ggsave(
  file.path(output_dir, paste0("cerno_treeplot_", tolower(plot_ontology), ".png")),
  p_tree,
  width = 12,
  height = 9,
  dpi = 300,
  limitsize = FALSE
)

# GO hierarchy graph ----------------------------------------------------------
p_go_graph <- cerno_go_graph(
  res_df = go_results,
  go_term_map = go_term_map,
  ontology = plot_ontology,
  top_n = top_n,
  padj_cutoff = padj_cutoff,
  selection_mode = "significant_or_top_n",
  max_depth = 4,
  include_roots = TRUE,
  label_ancestors = FALSE,
  layout = "sugiyama"
)

ggsave(
  file.path(output_dir, paste0("cerno_go_graph_", tolower(plot_ontology), ".png")),
  p_go_graph,
  width = 12,
  height = 9,
  dpi = 300,
  limitsize = FALSE
)

# GO semantic redundancy analysis and network --------------------------------
n_significant <- sum(
  go_results$ONTOLOGY == plot_ontology &
    !is.na(go_results$adj.P.Val) &
    go_results$adj.P.Val < padj_cutoff
)

redundancy_selection_mode <- if (n_significant > 0) {
  "significant_top_n"
} else {
  "top_n"
}

simplified_go <- cerno_simplify_go_terms(
  res_df = go_results,
  go_term_map = go_term_map,
  ontology = plot_ontology,
  padj_cutoff = padj_cutoff,
  top_n = top_n,
  selection_mode = redundancy_selection_mode,
  similarity_cutoff = 0.7,
  semantic_measure = "Wang"
)

write_csv(
  simplified_go$representatives,
  file.path(output_dir, paste0("GO_", plot_ontology, "_representatives.csv"))
)

write_csv(
  simplified_go$removed,
  file.path(output_dir, paste0("GO_", plot_ontology, "_redundant_terms.csv"))
)

p_redundancy <- cerno_go_redundancy_emap(
  simplify_result = simplified_go,
  label = "representatives",
  similarity_cutoff = 0.7,
  layout = "fr",
  show_removed_labels = TRUE
)

ggsave(
  file.path(
    output_dir,
    paste0("cerno_go_redundancy_", tolower(plot_ontology), ".png")
  ),
  p_redundancy,
  width = 12,
  height = 9,
  dpi = 300,
  limitsize = FALSE
)

# Optional Reactome plots -----------------------------------------------------
if (!is.null(reactome_results)) {
  p_reactome_dot <- cerno_dotplot(
    df = reactome_results,
    top_n = top_n,
    padj_cutoff = padj_cutoff,
    highlight1 = highlight_group_1,
    highlight2 = highlight_group_2,
    point_size_range = c(1.5, 7)
  )

  ggsave(
    file.path(output_dir, "cerno_dotplot_reactome.png"),
    p_reactome_dot,
    width = 10,
    height = 8,
    dpi = 300
  )

  p_reactome_tree <- cerno_treeplot(
    res_df = reactome_results,
    members_df = reactome_members,
    ontology = "REACTOME",
    top_n = top_n,
    padj_cutoff = padj_cutoff,
    selection_mode = "significant_or_top_n",
    nCluster = 5,
    show_cluster_labels = FALSE,
    show_legend = TRUE,
    split_after_words = 20
  )

  ggsave(
    file.path(output_dir, "cerno_treeplot_reactome.png"),
    p_reactome_tree,
    width = 12,
    height = 9,
    dpi = 300,
    limitsize = FALSE
  )
}
