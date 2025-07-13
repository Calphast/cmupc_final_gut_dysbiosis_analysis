# Gut Dysbiosis and Epithelial scRNA-seq Analysis

This repository contains the full R-based pipeline for my CMU Pre-College Computational Biology final project. The project analyzes epithelial cell responses to gut dysbiosis using two publicly available single-cell RNA sequencing (scRNA-seq) datasets from the Gene Expression Omnibus (GEO): **GSE168077** and **GSE169749**.

The main goals of this project are:

- To identify and annotate epithelial cell clusters
- To compute condition-specific cluster proportions and expansion ratios
- To explore how gut dysbiosis affects cellular composition and gene expression

---

## 📂 Data

This analysis uses two GEO datasets:

- **[GSE168077](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE168077)**: Gut epithelial single-cell transcriptomic profiles under dysbiosis
- **[GSE169749](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE169749)**: Keratinocyte/epithelial profiles from related conditions

**Note:** Due to file size, raw data is **not included** in this repository. To reproduce the analysis:

1. Visit the GEO links above and download the data manually
2. Place the extracted files in `data/GSE168077_RAW/` and `data/GSE169749_RAW/`

See `data/README.md` for more detailed instructions.

---

## ⚙️ How to Run the Analysis

### 1. Clone the Repository

```bash
git clone https://github.com/YOUR_USERNAME/gut-dysbiosis-analysis.git
cd gut-dysbiosis-analysis
```

### 2. Open the Project in RStudio

Open the file `gut_dysbiosis_analysis.Rproj`.

### 3. Install Required R Packages

You can install the necessary packages with:

```r
install.packages(c("Seurat", "tidyverse", "patchwork"))
```

If additional packages are needed (e.g. Bioconductor packages), use:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")
```

### 4. Run the Analysis Scripts

Scripts are organized into the `analysis/` folder:

- `analysis.R`: Primary analysis of GSE168077_RAW data
- `analysis2.R`: Primary analysis of GSE169749_RAW data (used in presentation)
- `keratinocyte_analysis.R`: Secondary analysis of GSE169749_RAW data focused on keratnocyte cells

---

## 📊 Results Summary

- UMAP plots showing cluster structure
- Proportion tables comparing conditions
- Cluster expansion ratios

---

## 🚫 Files Not Included

To keep the repository clean and lightweight, the following are **excluded**:

- Raw `.rds`, `.csv`, `.h5`, or `.gz` data files
- `.RData`, `.Rhistory`, and temporary workspace files
- Intermediate or large outputs

These are excluded via `.gitignore`.

---

## 👤 Author

**Gabriel Souto**  
CMU Pre-College Computational Biology  
Summer 2025
