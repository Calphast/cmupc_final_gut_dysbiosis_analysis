# Gut Dysbiosis and Epithelial scRNA-seq Analysis

This repository contains the full R-based pipeline for Pod23's CMU Pre-College Computational Biology final project. The project analyzes epithelial cell responses to gut dysbiosis using two publicly available single-cell RNA sequencing (scRNA-seq) datasets from the Gene Expression Omnibus (GEO): **GSE168077** and **GSE169749**.

The main goals of this project are:

- To identify and annotate cell clusters
- To compute condition-specific cluster proportions and expansion ratios
- To explore how gut dysbiosis affects cellular composition and gene expression

## 📂 Data

This analysis uses two GEO datasets:

- **[GSE168077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168077)**: Gut epithelial single-cell transcriptomic profiles under dysbiosis
- **[GSE169749](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169749)**: Keratinocyte/epithelial profiles from related conditions

**Note:** Due to file size, raw data is **not included** in this repository. To reproduce the analysis:

1. Visit the GEO links above and download the data manually
2. Place the extracted files in `data/GSE168077_RAW/` and `data/GSE169749_RAW/`

Scripts are organized into the `analysis/` folder:

- `analysis.R`: Primary analysis of GSE168077_RAW data
- `analysis2.R`: Primary analysis of GSE169749_RAW data (used in presentation)
- `keratinocyte_analysis.R`: Secondary analysis of GSE169749_RAW data focused on keratnocyte cells

## 🔧 Dependencies

This project was developed using the following R packages:

- `Seurat` (v5 or compatible)
- `tidyverse`
- `patchwork`
- `ggplot2`
- `dplyr`

## 📊 Results Summary

- UMAP plots showing cluster structure
- Proportion tables comparing conditions
- Cluster expansion ratios

## 🚫 Files Not Included

To keep the repository clean and lightweight, the following are **excluded**:

- Raw `.rds`, `.csv`, `.h5`, or `.gz` data files
- `.RData`, `.Rhistory`, and temporary workspace files
- Intermediate or large outputs

These are excluded via `.gitignore`.

## 👤 Author

**Gabriel Souto**  
CMU Pre-College Computational Biology  
Summer 2025
