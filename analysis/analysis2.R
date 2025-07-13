library(Seurat)
library(ggplot2)
library(dplyr)

gut_healthy <- Load10X_Spatial(
  data.dir = "data/GSE169749_RAW/A1",
  filename = "filtered_feature_bc_matrix.h5",
  image = "tissue_hires_image.png")

gut_healthy$conditions <- "Healthy"

gut_dss <- Load10X_Spatial(
  data.dir = "data/GSE169749_RAW/B1",
  filename = "filtered_feature_bc_matrix.h5",
  image = "tissue_hires_image.png.gz"
)

gut_dss$conditions <- "DSS"

# cleaning up Healthy dataset
head(rownames(gut_healthy))

gut_healthy[["percent.mt"]] <- PercentageFeatureSet(gut_healthy, pattern = "^mt-")

DefaultAssay(gut_healthy) <- "Spatial"  # or "Spatial" if using NormalizeData()
VlnPlot(gut_healthy, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), slot = "counts")

gut_healthy <- subset(gut_healthy, subset = nFeature_Spatial > 1250 & nCount_Spatial < 41000 & percent.mt < 10)

gut_healthy <- SCTransform(gut_healthy,
                     assay = "Spatial",
                     return.only.var.genes = FALSE,
                     verbose = FALSE)

gut_healthy <- RunPCA(gut_healthy, assay = "SCT", verbose = FALSE)
gut_healthy <- FindNeighbors(gut_healthy, dims = 1:30)
gut_healthy <- FindClusters(gut_healthy, resolution = 0.4)
gut_healthy <- RunUMAP(gut_healthy, dims = 1:30)

DimPlot(gut_healthy, reduction = "umap", label = TRUE)

# cleaning up DSS dataset
gut_dss[["percent.mt"]] <- PercentageFeatureSet(gut_dss, pattern = "^mt-")

DefaultAssay(gut_dss) <- "Spatial"  # or "Spatial" if using NormalizeData()
VlnPlot(gut_dss, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), slot = "counts")

gut_dss <- subset(gut_dss, subset = nFeature_Spatial > 2000 & nCount_Spatial < 30000 & percent.mt < 9)

gut_dss <- SCTransform(gut_dss,
                  assay = "Spatial",
                  return.only.var.genes = FALSE,
                  verbose = FALSE)

gut_dss <- RunPCA(gut_dss, assay = "SCT", verbose = FALSE)
gut_dss <- FindNeighbors(gut_dss, dims = 1:30)
gut_dss <- FindClusters(gut_dss, resolution = 0.4)
gut_dss <- RunUMAP(gut_dss, dims = 1:30)

DimPlot(gut_dss, reduction = "umap", label = TRUE)

# 3. prep & integrate using the same list
DefaultAssay(gut_healthy) <- "SCT"
DefaultAssay(gut_dss)     <- "SCT"

gut_list <- list(gut_healthy, gut_dss)
features <- SelectIntegrationFeatures(
  object.list = gut_list,
  nfeatures   = 3000
)

gut_list <- PrepSCTIntegration(object.list = gut_list, anchor.features = features)

anchors <- FindIntegrationAnchors(
  object.list          = gut_list,      # use the list, not re-specify two objects
  normalization.method = "SCT",
  anchor.features      = features
)

gut_combined <- IntegrateData(
  anchorset            = anchors,
  normalization.method = "SCT"
)

# Integration Analysis
gut_combined <- RunPCA(gut_combined, verbose = FALSE)

gut_combined <- FindNeighbors(gut_combined, dims = 1:30)
gut_combined <- FindClusters(gut_combined, resolution = 0.4)

gut_combined <- RunUMAP(gut_combined, dims = 1:30)

DimPlot(gut_combined, reduction = "pca",   label = TRUE)    # PCA
DimPlot(gut_combined, reduction = "umap",  label = TRUE)    # UMAP

# Getting the markers of each cell
markers <- FindAllMarkers(
  object          = gut_combined,
  only.pos        = TRUE,
  min.pct         = 0.25,
  logfc.threshold = 0.25
)

top_markers <- markers %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = pct.1 - pct.2) %>%
  group_by(cluster) %>%
  slice_max(order_by = pct_diff, n = 5)

new.ids <- c(
  "0" = "Endothelial Cells",
  "1" = "Fibroblasts",
  "2" = "Enterocytes",
  "3" = "Neurons",
  "4" = "B Cells",
  "5" = "Enterocytes",
  "6" = "Enterocytes",
  "7" = "Neutrophils",
  "8" = "Keratinocytes",
  "9" = "Neurons",
  "10" = "Macrophages"
)

gut_combined <- RenameIdents(gut_combined, new.ids)
gut_combined$seurat_annotations <- Idents(gut_combined)

DimPlot(gut_combined, reduction = "umap", group.by = c("stim", "seurat_annotations")) # Visualizing stim and celltype

# additional analysis to see what is expanding and what is contracting
# within the cluster of normal vs compromised
table(Idents(gut_combined), gut_combined$conditions)
prop.table(table(Idents(gut_combined), gut_combined$conditions), margin = 1)

ct <- table(
  Cluster   = Idents(gut_combined),
  Condition = gut_combined$conditions
)

prop_cluster <- prop.table(ct, margin = 1)

expansion_ratio <- prop_cluster[ , "DSS"] / prop_cluster[ , "Healthy"]

expansion_ratio

# Create contingency table
ct <- table(Idents(gut_combined), gut_combined$conditions)

# Compute proportion table (row-wise)
proportions <- prop.table(ct, margin = 1)

# Compute expansion ratio per cluster (treated / control)
expansion_ratio <- ct[, "Healthy"] / ct[, "DSS"]

# Combine all into one data frame
result <- cbind(
  Control_Count = ct[, "Healthy"],
  Treated_Count = ct[, "DSS"],
  Treated_Proportion = round(proportions[, "DSS"], 3),
  Control_Proportion = round(proportions[, "Healthy"], 3),
  Expansion_Ratio = round(expansion_ratio, 2)
)

# Convert to data frame for plotting or presentation
result_df <- as.data.frame(result)

result_df
