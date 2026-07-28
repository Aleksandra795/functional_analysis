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

## Key parameters

The visualization functions use a consistent set of parameters wherever applicable:

* `highlight1` and `highlight2` are lists of exact pathway names or IDs that will be highlighted in the dotplot (by color) and treeplot (with one or two stars *, respectively), optional.
* `ontology` selects the functional category to analyze. Supported GO values are `"BP"`, `"CC"` and `"MF"`. The tree plot also supports `"REACTOME"` when appropriate pathway definitions are supplied.
* `top_n` defines the maximum number of terms included in the visualization.
* `padj_cutoff` defines the adjusted p-value threshold used to identify statistically significant terms. The default value is `0.05`.
* `selection_mode` controls how terms are selected before constructing tree, hierarchy and redundancy visualizations:

  * `"significant_top_n"` selects up to `top_n` terms with adjusted p-values below `padj_cutoff`.
  * `"significant_all"` selects all terms with adjusted p-values below `padj_cutoff`.
  * `"top_n"` selects the `top_n` terms with the lowest adjusted p-values, regardless of significance.
  * `"significant_or_top_n"` uses significant terms when enough are available for the selected visualization; otherwise, it falls back to the `top_n` terms with the lowest adjusted p-values.
* `ontology_col`, `id_col`, `term_col`, `padj_col` and `auc_col` specify the corresponding columns in the input data. Their defaults are `"ONTOLOGY"`, `"ID"`, `"Title"`, `"adj.P.Val"` and `"AUC"`.
* `split_after_words` controls line wrapping in long pathway names.
* `point_size_range` controls the minimum and maximum point sizes used to represent AUC values.
* `nCluster` defines the number of functional clusters displayed in a pathway tree.
* `similarity_cutoff` defines the minimum GO semantic similarity required to group or connect terms in the redundancy analysis.
* `semantic_measure` selects the GO semantic similarity method. The default is `"Wang"`.
* `max_depth` controls how many GO ancestor levels are included in the hierarchy graph. For greater clarity, if you have >30 significant pathways, set max_depth to 2-3.

And other parameters related to labels, legends, layout and figure dimensions.

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
