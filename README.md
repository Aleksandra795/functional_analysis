# functional_analysis

R functions for visualizing functional enrichment results obtained with the CERNO test from tmod package. The repository supports Gene Ontology Biological Process, Cellular Component and Molecular Function collections, as well as Reactome pathways.

`CERNO_vis.R` is the main script. Set the input directory, output directory and input file names in its **User settings** section, then run the script from the repository root.

## Installation

R 4.1 or later is recommended.

```r
source("install_packages.R")
```

## Required input columns

CERNO results may be stored in Excel files and should contain the following columns:

| Column | Meaning | Required by |
|---|---|---|
| `ID` | Pathway or gene-set identifier | tree plot, GO graph, redundancy analysis |
| `Title` | Pathway name | all visualizations |
| `adj.P.Val` | FDR-adjusted p-value | all visualizations |
| `AUC` | CERNO effect-size measure | all visualizations |
| `ONTOLOGY` | `BP`, `CC`, `MF` or `REACTOME` | added automatically by `CERNO_vis.R` |
| `n_pathway_total` | Pathway size | optional; used by filtering helpers |

`ID` should preferably correspond to an MSigDB gene-set identifier. The main script normalizes IDs against the collections returned by `msigdbr`.

## Main functions

| Function | Purpose | Most important parameters |
|---|---|---|
| `cerno_dotplot()` | Displays pathways by `-log10(FDR)` with point size representing AUC | `top_n`, `padj_cutoff`, `highlight1`, `highlight2`, `point_size_range` |
| `cerno_manhattan_go()` | Compares GO terms across BP, CC and MF | `padj_cutoff`, `label_top_n`, `label_padj_cutoff`, `max_terms_per_ontology`, `jitter_width` |
| `cerno_treeplot()` | Hierarchically clusters pathways using Jaccard overlap of member genes | `ontology`, `top_n`, `selection_mode`, `nCluster`, `cluster_method` |
| `cerno_go_graph()` | Displays selected GO terms and their ancestors in the GO hierarchy | `ontology`, `top_n`, `selection_mode`, `max_depth`, `include_roots`, `layout` |
| `cerno_simplify_go_terms()` | Groups semantically redundant GO terms and selects representatives | `ontology`, `selection_mode`, `similarity_cutoff`, `semantic_measure` |
| `cerno_go_redundancy_emap()` | Visualizes semantic similarity among GO terms | `label`, `similarity_cutoff`, `layout`, `show_removed_labels` |

Common values of `selection_mode` are `significant_top_n`, `significant_all`, `top_n` and `significant_or_top_n`, depending on the function.

## Running the analysis

```r
source("CERNO_vis.R")
```

The script creates PNG figures, CSV files describing representative and redundant GO terms in the selected output directory.

## Acknowledgements

The visualization concepts and analytical workflows implemented in this
repository were inspired by the functional enrichment visualization methods
presented in the Biomedical Knowledge Mining book and the enrichplot
ecosystem.

The code provides custom implementations adapted to tabular CERNO results.
In particular, the dot plot, Manhattan-style plot, pathway-overlap tree,
GO hierarchy graph and semantic redundancy network were adapted to operate
on CERNO result tables.

CERNO results may be generated using the tmod R package.

References:

- Biomedical Knowledge Mining:
  https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
- enrichplot:
  https://bioconductor.org/packages/enrichplot/
- tmod:
  https://CRAN.R-project.org/package=tmod

This repository is not affiliated with or endorsed by the authors of
enrichplot, clusterProfiler or tmod.

## License

This project is distributed under the GNU General Public License v3.0.
See the `LICENSE` file for details.
