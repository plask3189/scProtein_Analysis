# scProtein_Analysis

An R pipeline for single-cell multi-omic analysis of DNA variants and surface protein expression from Mission Bio Tapestri data. Originally developed to study **BRAF-mutant AML** samples, with support for cohort-level comparisons.

## Overview

This pipeline processes `.h5` files from the Tapestri platform and performs:

1. **Variant identification** — Calls somatic variants (e.g. FLT3-ITD, NRAS, KRAS, PTPN11, KIT, DNMT3A, TET2, IDH1/2, NPM1, SRSF2, ASXL1) per cell using configurable genotype quality and VAF thresholds.
2. **Protein normalization** — Applies DSB (Denoised and Scaled by Background) and CLR normalization to antibody-derived tag (ADT) counts, using empty droplets as a background reference.
3. **Clustering** — Builds a Seurat object from normalized protein data and performs dimensionality reduction and clustering.
4. **Visualization** — Generates UMAPs, violin plots, heatmaps, and significance bracket figures.
5. **Cohort comparison** (optional) — Merges samples from two groups (` | **Main entry point.** Orchestrates the full pipeline from `.h5` input through normalization, clustering, and saving results. |
| `braf_protein_hng-marker handling. |
| `load_libraries.R` | C(PCA, UMAP, Louvain). |
| `protein_normaion.R` | Differential protein expression teerring variant/clone labels onto Seurat objects. |
| `violinor genotypes. |
| `doheat
