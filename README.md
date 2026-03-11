# scProtein_Analysis

R code for single-cell multiomic analysis of DNA variants and surface protein expression from Mission Bio Tapestri data. Generates Figure 2: Genotype-immunophenotype relationships in BRAF-mutant AML clones and Supp. Figure 1.  Uses functions from: 
https://github.com/bowmanr/scDNA 
scDNA: Single Cell DNA analysis software toolkit for subclonality discovery and assessment Michael Bowman, Shreeya Gounder, Varsha Singh, Olga Shestova, Troy Robinson, Amy Zhang, Anushka Gandhi, Roopsha Bandopadhyay, Sheng F. Cai, Ross L. Levine, Saar I. Gill, Linde A. Miles, Robert L. Bowman bioRxiv 2025.12.19.694255; doi: https://doi.org/10.64898/2025.12.19.694255


## Overview

Processes `.h5` files from the Tapestri platform and performs:

1. **Variant identification** — Calls somatic variants (e.g. FLT3-ITD, NRAS, KRAS, PTPN11, KIT, DNMT3A, TET2, IDH1/2, NPM1, SRSF2, ASXL1) per cell using configurable genotype quality and VAF thresholds.
2. **Protein normalization** — Applies DSB (Denoised and Scaled by Background) and CLR normalization to antibody-derived tag (ADT) counts, using empty droplets as a background reference.
3. **Clustering** — Builds a Seurat object from normalized protein data and variant metadata. Performs dimensionality reduction and clustering
4. **Visualization** — Generates UMAPs, violin plots, heatmaps, and significance analysis.
