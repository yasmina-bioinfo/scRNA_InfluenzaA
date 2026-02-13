# 10_Tcells_Subset.R
# Subset strict T cells (exclude NK) using:
# cluster-based prefilter + CD3 gene expression safety filter
# Then re-run standard Seurat workflow on the subset and compute markers.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# 0) Load object (Harmony merged)
# -----------------------------
seurat_rds <- "Results/seurat_objects/seu_harmony_merged.rds"
if (!file.exists(seurat_rds)) stop("Seurat RDS not found: ", seurat_rds)

seu <- readRDS(seurat_rds)

# -----------------------------
# 1) Parameters you may EDIT
# -----------------------------
# Cluster-based prefilter: put here the global clusters that contain T cells
# (start with c("3","5") based on your current annotation; adjust if needed)
tcell_clusters_prefilter <- c("3", "5")

# CD3 safety genes (at least one must be detected per cell)
cd3_genes <- c("CD3D", "CD3E", "TRAC")

# Reclustering params for T cell subset
nfeatures_hvg <- 2000
dims_use <- 1:30
resolution_use <- 0.6

# Output paths
out_dir <- "Results/tcells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 2) Safety checks
# -----------------------------
if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("Column 'seurat_clusters' not found in meta.data.")
}
seu$seurat_clusters <- as.factor(seu$seurat_clusters)

# Check CD3 genes exist
missing_cd3 <- cd3_genes[!cd3_genes %in% rownames(seu)]
if (length(missing_cd3) == length(cd3_genes)) {
  stop("None of the CD3 genes were found in the object: ", paste(cd3_genes, collapse = ", "))
}

# -----------------------------
# 3) Step C: cluster-based prefilter
# -----------------------------
prefilter_cells <- WhichCells(seu, idents = tcell_clusters_prefilter)

# If Idents() isn't set to seurat_clusters, set it temporarily
Idents(seu) <- seu$seurat_clusters
prefilter_cells <- WhichCells(seu, idents = tcell_clusters_prefilter)

seu_pref <- subset(seu, cells = prefilter_cells)

# -----------------------------
# 4) Step C: CD3 expression safety filter
# Keep cells with at least one of CD3D/CD3E/TRAC detected (>0)
# -----------------------------
present_cd3 <- cd3_genes[cd3_genes %in% rownames(seu_pref)]

present_cd3 <- cd3_genes[cd3_genes %in% rownames(seu_pref)]
if (length(present_cd3) == 0) stop("No CD3 genes found in subset object.")

# FetchData returns a data.frame (cells x genes) using the default assay/layer
cd3_df <- FetchData(seu_pref, vars = present_cd3)

# Keep cells where at least one CD3 gene is detected (>0)
keep_cells <- rownames(cd3_df)[rowSums(cd3_df > 0) >= 1]

seu_t <- subset(seu_pref, cells = keep_cells)

# Report counts
message("Prefilter (clusters ", paste(tcell_clusters_prefilter, collapse = ", "), "): ", ncol(seu_pref), " cells")
message("After CD3 safety filter: ", ncol(seu_t), " cells")

# Save a small table of filtering steps
filter_tbl <- data.frame(
  step = c("prefilter_clusters", "cd3_safety_filter"),
  n_cells = c(ncol(seu_pref), ncol(seu_t))
)
write.table(filter_tbl, file = file.path(out_dir, "tcell_filtering_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# -----------------------------
# 5) Re-run Seurat workflow on T cell subset
# -----------------------------
DefaultAssay(seu_t) <- "RNA"
seu_t <- NormalizeData(seu_t)
seu_t <- FindVariableFeatures(seu_t, selection.method = "vst", nfeatures = nfeatures_hvg)
seu_t <- ScaleData(seu_t, features = VariableFeatures(seu_t))
seu_t <- RunPCA(seu_t, features = VariableFeatures(seu_t))
seu_t <- FindNeighbors(seu_t, dims = dims_use)
seu_t <- FindClusters(seu_t, resolution = resolution_use)
seu_t <- RunUMAP(seu_t, dims = dims_use)

# -----------------------------
# 6) Markers on T cell subset
# -----------------------------

# --- Seurat v5: Join layers (needed for DE / FindAllMarkers) ---
# Do this on the T-cell subset only (smaller -> RAM OK)
if ("RNA" %in% names(seu_t@assays)) {
  # Join layers inside the RNA assay
  seu_t[["RNA"]] <- JoinLayers(seu_t[["RNA"]])
} else {
  stop("RNA assay not found in seu_t.")
}

markers_t <- FindAllMarkers(
  seu_t,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.table(markers_t, file = file.path(out_dir, "markers_tcells_all_clusters.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Top20 per T-cell cluster (condensed) - robust to column names
cluster_col <- if ("cluster" %in% colnames(markers_t)) {
  "cluster"
} else if ("ident" %in% colnames(markers_t)) {
  "ident"
} else {
  stop("No cluster/ident column found in markers_t. Columns are: ",
       paste(colnames(markers_t), collapse = ", "))
}

gene_col <- if ("gene" %in% colnames(markers_t)) "gene" else NULL
if (is.null(gene_col)) stop("No 'gene' column found in markers_t.")

summary_t <- markers_t %>%
  group_by(.data[[cluster_col]]) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  summarise(top_genes = paste(.data[[gene_col]], collapse = ", "), .groups = "drop") %>%
  rename(cluster = .data[[cluster_col]]) %>%
  arrange(as.numeric(as.character(cluster)))

# -----------------------------
# 7) UMAP plot for T cell subset
# -----------------------------
p_t_umap <- DimPlot(seu_t, reduction = "umap", label = TRUE, repel = TRUE) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_blank())

ggplot2::ggsave(
  filename = file.path(out_dir, "UMAP_tcells_by_cluster.png"),
  plot = p_t_umap,
  width = 6, height = 6, units = "in", dpi = 600
)

# Save object
saveRDS(seu_t, file = file.path(out_dir, "seu_tcells_subset_reclustered.rds"))

message("\nDONE ✅ Saved in: ", out_dir, "\n",
        "- tcell_filtering_summary.tsv\n",
        "- markers_tcells_all_clusters.tsv\n",
        "- top20_markers_per_tcell_cluster.tsv\n",
        "- UMAP_tcells_by_cluster.png\n",
        "- seu_tcells_subset_reclustered.rds\n")
