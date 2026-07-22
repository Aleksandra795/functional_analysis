cran_packages <- c(
  "ape",
  "dplyr",
  "ggraph",
  "ggplot2",
  "ggrepel",
  "ggtext",
  "igraph",
  "msigdbr",
  "purrr",
  "readr",
  "readxl",
  "stringr",
  "tibble"
)

missing_cran <- cran_packages[
  !vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_cran) > 0) {
  install.packages(missing_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c(
  "AnnotationDbi",
  "GO.db",
  "GOSemSim",
  "ggtree",
  "org.Hs.eg.db"
)

missing_bioc <- bioc_packages[
  !vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}
