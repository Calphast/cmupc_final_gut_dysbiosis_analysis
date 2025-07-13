source("analysis2.R")

ker_subset_strict <- subset(
  gut_combined,
  subset = Krt6b  > 0 &
    Tacstd2 > 0 &
    Calml3  > 0 &
    Prss22  > 0 &
    Arg1    > 0 & 
    seurat_annotations == "Keratinocytes"
)

ker_subset_strict

gut_combined <- AddModuleScore(gut_combined, 
                               features = list(c("Krt6b", "Tacstd2", "Calml3", "Prss22", "Arg1")),
                               name = "ker_score")

VlnPlot(gut_combined, features = "ker_score1", pt.size = 0.1)

ker_subset_score <- subset(gut_combined, subset = ker_score1 > 0.0)

VlnPlot(ker_subset_score, features = c("percent.mt","nFeature_RNA"))

ker_subset_score <- SCTransform(ker_subset_score, ) %>%
  RunPCA() %>%
  FindNeighbors(dims = 1:20) %>%
  FindClusters() %>%
  RunUMAP(dims = 1:20)
