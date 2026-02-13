# Scripts/08_Markers_Global_downsample.R
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

seu <- readRDS("Results/seurat_objects/seu_harmony_merged.rds")
DefaultAssay(seu) <- "RNA"
Idents(seu) <- "seurat_clusters"

# ---- Downsample: max N cells per cluster (adjustable) ----
max_cells_per_cluster <- 1500  # 1000–2000 is usually safe on laptop

set.seed(1)
cells_keep <- unlist(lapply(levels(Idents(seu)), function(cl) {
  cells <- WhichCells(seu, idents = cl)
  if (length(cells) > max_cells_per_cluster) sample(cells, max_cells_per_cluster) else cells
}))

seu_ds <- subset(seu, cells = cells_keep)
rm(seu); gc(FALSE)

# Join layers on the downsampled object (now should fit in RAM)
seu_ds[["RNA"]] <- JoinLayers(seu_ds[["RNA"]])

# Markers
dir.create("Results/markers", recursive = TRUE, showWarnings = FALSE)

markers <- FindAllMarkers(
  seu_ds,
  group.by = "seurat_clusters",
  assay = "RNA",
  slot = "data",
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.10,
  return.thresh = 1
)

# Save full table
write.table(markers, "Results/markers/markers_all_clusters_downsample.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Print top 20 per cluster
top_n <- 20
top_markers <- markers %>%
  mutate(gene = ifelse("gene" %in% colnames(markers), gene, rownames(markers))) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = top_n) %>%
  select(cluster, gene, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
  arrange(as.numeric(as.character(cluster)), desc(avg_log2FC))

message("\n=== TOP ", top_n, " MARKERS PER CLUSTER (downsampled) ===\n")
print(top_markers, n = Inf)

message("\nDONE: downsample markers saved -> Results/markers/markers_all_clusters_downsample.tsv")

# ---- Condensed summary (one line per cluster) ----
summary_clusters <- markers %>%
  mutate(gene = ifelse("gene" %in% colnames(markers), gene, rownames(markers))) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop") %>%
  arrange(as.numeric(as.character(cluster)))

write.table(summary_clusters,
            "Results/annotation/top20_markers_per_cluster.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\n=== SUMMARY (Top 20 genes per cluster) ===\n")
print(summary_clusters, n = Inf)

