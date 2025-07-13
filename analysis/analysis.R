library(Seurat)
library(SeuratData)
library(patchwork)

# instantiating the data and merging it into one integrated object
wt.data <- Read10X(data.dir = "data/GSE168077_RAW/WT")
wt <- CreateSeuratObject(counts = wt.data, project = "WT1")
wt$stim <- "WT"
cko.data <- Read10X(data.dir = "data/GSE168077_RAW/cKO")
cko <- CreateSeuratObject(counts = cko.data, project = "cKO1")
cko$stim <- "cKO"

gda <- merge(wt, y = cko, add.cell.ids = c("WT", "CKO"), project = "GDA1")

# Analysis without integration
# run standard anlaysis workflow
gda <- NormalizeData(gda)
gda <- FindVariableFeatures(gda)
gda <- ScaleData(gda)
gda <- RunPCA(gda)

gda <- FindNeighbors(gda, dims = 1:30, reduction = "pca")
gda <- FindClusters(gda, resolution = 2, cluster.name = "unintegrated_clusters")

gda <- RunUMAP(gda, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(gda, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

# performing integration
gda <- IntegrateLayers(object = gda, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                       verbose = FALSE)

# re-join layers after integration
gda[["RNA"]] <- JoinLayers(gda[["RNA"]])

gda <- FindNeighbors(gda, reduction = "integrated.cca", dims = 1:30)
gda <- FindClusters(gda, resolution = 1)

gda <- RunUMAP(gda, dims = 1:30, reduction = "integrated.cca")
DimPlot(gda, reduction = "umap", group.by = c("stim", "seurat_annotations"))

# Getting the markers of each cell
markers <- FindAllMarkers(
  object          = gda,
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
  "0" = "NK Cells",
  "1" = "Neurons", 
  "2" = "Neurons",
  "3"= "Macrophages",
  "4" = "Keratinocytes",
  "5" = "Keratinocytes",
  "6" = "Endothelial Cells",
  "7" = "Endothelial Cells", 
  "8" = "Nuocytes", 
  "9" = "Fibroblasts",
  "10" = "Fibroblasts", 
  "11" = "Endothelial Cells",
  "12" = "Endothelial Cells", 
  "13" = "T cells",
  "14" = "NK cells",
  "15" = "Endothelial Cells",
  "16" = "T cells",
  "17" = "Endothelial Cells",
  "18" = "T cells"
)

gda <- RenameIdents(gda, new.ids)
gda$seurat_annotations <- Idents(gda)

DimPlot(gda, reduction = "umap", group.by = c("stim", "seurat_annotations")) # Visualizing stim and celltype
DimPlot(gda, reduction = "umap", split.by = "stim", group.by = "seurat_annotations") # Visualizing conditions side-by-side
