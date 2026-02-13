# 07_Figures: PCA and UMAP visualization after Harmony batch correction

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# Load Harmony object
seu <- readRDS("Results/seurat_objects/seu_harmony_merged.rds")

# Create output directory
outdir <- "Results/figures/harmony"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---- PCA ----
p_pca_clusters <- DimPlot(
  seu,
  reduction = "pca",
  group.by = "seurat_clusters",
  raster = TRUE
) + ggtitle("PCA - Seurat clusters")

p_pca_samples <- DimPlot(
  seu,
  reduction = "pca",
  group.by = "sample",
  raster = TRUE
) + ggtitle("PCA - Samples")

# ---- UMAP ----
p_umap_clusters <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = TRUE
) + ggtitle("UMAP - Seurat clusters")

p_umap_samples <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "sample",
  raster = TRUE
) + ggtitle("UMAP - Samples")

# ---- Save figures ----
ggsave(file.path(outdir, "PCA_clusters.png"),  p_pca_clusters, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "PCA_samples.png"),   p_pca_samples,  width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "UMAP_clusters.png"), p_umap_clusters, width = 8, height = 6, dpi = 300)
ggsave(file.path(outdir, "UMAP_samples.png"),  p_umap_samples,  width = 8, height = 6, dpi = 300)

message("DONE: PCA and UMAP figures saved.")
