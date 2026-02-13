# 09_Celltype_Annotation.R
# Purpose: add broad Cell_Annotation labels (family-level) from seurat_clusters
#          and generate a global UMAP colored by Cell_Annotation.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# -----------------------------
# 0) Load Seurat object
# -----------------------------
# EDIT THIS: path to your integrated object (Harmony) saved as .rds
seurat_rds <- "Results/seurat_objects/seu_harmony_merged.rds"

if (!file.exists(seurat_rds)) {
  stop("Seurat RDS not found: ", seurat_rds)
}

seu <- readRDS(seurat_rds)

# -----------------------------
# 1) Safety checks
# -----------------------------
if (!"seurat_clusters" %in% colnames(seu@meta.data)) {
  stop("Column 'seurat_clusters' not found in meta.data.")
}

# Ensure it's factor (you already checked, but keep it robust)
seu@meta.data$seurat_clusters <- as.factor(seu@meta.data$seurat_clusters)

# -----------------------------
# 2) Define broad family mapping (0–22)
# -----------------------------
# Based on your validated labels:
# Myeloid: 0,2,7,8,17,19,20,21
# Neutrophils: 1,14
# Lymphoid: 3,4,5,9,11,13
# Mast/Baso: 16,22
# Platelets: 12
# RBC: 10
# Cycling: 18
# Unknown/Rare non-immune: 6,15
cluster_to_family <- c(
  "0"  = "Monocytes/Macrophages/DC",
  "1"  = "Neutrophiles",
  "2"  = "Monocytes/Macrophages/DC",
  "3"  = "Lymphocytes T/NK",
  "4"  = "Lymphocytes T/NK",
  "5"  = "Lymphocytes T/NK",
  "6"  = "Others",
  "7"  = "Monocytes/Macrophages/DC",
  "8"  = "Monocytes/Macrophages/DC",
  "9"  = "Lymphocytes B/Plasmocytes",
  "10" = "Others",
  "11" = "Lymphocytes B/Plasmocytes",
  "12" = "Others",
  "13" = "Lymphocytes T/NK",
  "14" = "Neutrophiles",
  "15" = "Others",
  "16" = "Others",
  "17" = "Monocytes/Macrophages/DC",
  "18" = "Others",   # Cycling integrated in "others" for this view
  "19" = "Monocytes/Macrophages/DC",
  "20" = "Monocytes/Macrophages/DC",
  "21" = "Monocytes/Macrophages/DC",
  "22" = "Others"
)

custom_colors <- c(
  "Lymphocytes T/NK"            = "#1F77B4",  # blue
  "Lymphocytes B/Plasmocytes"   = "#FFD700",  # yellow
  "Monocytes/Macrophages/DC"    = "#2CA02C",  # green
  "Neutrophiles"                = "#D62728",  # red
  "Others"                      = "#7F7F7F"   # grey
)

# Add Cell_Annotation
seu$Cell_Annotation <- unname(cluster_to_family[as.character(seu$seurat_clusters)])

# -----------------------------
# 3) Validate: no missing annotations
# -----------------------------
missing_clusters <- sort(unique(as.character(seu$seurat_clusters[is.na(seu$Cell_Annotation)])))
if (length(missing_clusters) > 0) {
  stop("Missing Cell_Annotation for clusters: ", paste(missing_clusters, collapse = ", "),
       "\nFix cluster_to_family mapping.")
}

# Quick summary table (saved + printed)
dir.create("Results/annotation", recursive = TRUE, showWarnings = FALSE)

summary_tbl <- seu@meta.data %>%
  count(seurat_clusters, Cell_Annotation, name = "n_cells") %>%
  arrange(as.numeric(as.character(seurat_clusters)))

write.table(summary_tbl,
            file = "Results/annotation/cell_annotation_summary.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\n=== Cell_Annotation summary ===\n")
print(as.data.frame(summary_tbl))

# -----------------------------
# 4) Global UMAP colored by Cell_Annotation
# -----------------------------
# Make sure UMAP exists
if (!"umap" %in% names(seu@reductions)) {
  stop("UMAP reduction not found in this object. Run UMAP before this script.")
}

p1 <- DimPlot(
  seu,
  reduction = "umap",
  group.by = "Cell_Annotation",
  cols = custom_colors,
  pt.size = 0.2
) +
  theme_classic(base_size = 11) +
  theme(
  legend.title = element_blank(),
  legend.text  = element_text(size = 12),
  legend.key.size = unit(0.6, "lines"),
  axis.title   = element_text(size = 9),
  axis.text    = element_text(size = 8),
  plot.title   = element_blank()
)

# Save plots
ggplot2::ggsave(
  "Results/annotation/UMAP_global_by_Cell_Annotation.pdf",
  plot = p1,
  width = 6,
  height = 6,
  units = "in"
)

ggplot2::ggsave(
  "Results/annotation/UMAP_global_by_Cell_Annotation.png",
  plot = p1,
  width = 6,
  height = 6,
  units = "in",
  dpi = 600
)

# Save updated object
saveRDS(seu, file = "Results/seurat/seurat_harmony_integrated_with_CellAnnotation.rds")

message("\nDONE ✅  Saved:\n",
        "- Results/annotation/cell_annotation_summary.tsv\n",
        "- Results/annotation/UMAP_global_by_Cell_Annotation.png\n",
        "- Results/annotation/UMAP_global_by_Cluster.png\n",
        "- Results/seurat/seurat_harmony_integrated_with_CellAnnotation.rds\n")
