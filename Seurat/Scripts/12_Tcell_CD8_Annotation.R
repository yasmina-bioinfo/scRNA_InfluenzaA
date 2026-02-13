#!/usr/bin/env Rscript

# ============================================================
# 12_Tcell_CD8_Annotated.R
# Goal:
#   - Load T-cell annotated object (Results/tcells/seu_tcells_annotated.rds)
#   - Subset CD8 populations
#   - Export avg expression table per CD8 cluster (0/4/6)
#   - Generate DotPlot (CD8 clusters)
#   - Rename clusters officially (CD8_Annotation)
#   - Generate UMAP (CD8_Annotation)
#   - Save CD8 annotated object
#
# Outputs:
#   - Results/tcells/CD8_cluster_avgexpr_markers.tsv
#   - Results/tcells/DotPlot_CD8_by_cluster.png
#   - Results/tcells/UMAP_CD8_annotated.png
#   - Results/tcells/seu_cd8_annotated.rds
# ============================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# ----------------------------
# Paths
# ----------------------------
in_rds  <- "Results/tcells/seu_tcells_annotated.rds"
out_dir <- "Results/tcells"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ----------------------------
# Load object
# ----------------------------
if (!file.exists(in_rds)) stop("Missing input RDS: ", in_rds)
seu_t <- readRDS(in_rds)

# ----------------------------
# Subset CD8 cells (based on your GLOBAL Tcell_Annotation)
# ----------------------------
if (!("Tcell_Annotation" %in% colnames(seu_t@meta.data))) {
  stop("Column 'Tcell_Annotation' not found in seu_t@meta.data. Did you save the annotated T-cell object?")
}

cd8_keep <- c("CD8_Effector", "CD8_Memory", "CD8_Terminal_Cytotoxic")
seu_cd8 <- subset(seu_t, subset = Tcell_Annotation %in% cd8_keep)

# Safety checks
if (ncol(seu_cd8) == 0) stop("CD8 subset is empty. Check cd8_keep labels match your Tcell_Annotation values.")
if (!("seurat_clusters" %in% colnames(seu_cd8@meta.data))) stop("seurat_clusters not found in CD8 object.")

cd8_clusters <- sort(unique(as.character(seu_cd8$seurat_clusters)))
message("CD8 clusters found: ", paste(cd8_clusters, collapse = ", "))

# ----------------------------
# Marker panel for CD8 dotplot/table
# (enough info, but not too many)
# ----------------------------
features_cd8 <- c(
  "CD8A","CD8B",
  "NKG7","GNLY","GZMB","PRF1","FGFBP2",
  "CCR7","IL7R",
  "PDCD1","LAG3","TIGIT","TOX",
  "IFNG","STAT1","ISG15",
  "MKI67"
)

features_present <- features_cd8[features_cd8 %in% rownames(seu_cd8)]
if (length(features_present) == 0) stop("None of the CD8 marker features are present in rownames(seu_cd8).")

# ----------------------------
# Avg expression per cluster table (CD8)
# Using AverageExpression (Seurat v5-safe)
# ----------------------------
Idents(seu_cd8) <- "seurat_clusters"

avg_list <- AverageExpression(
  seu_cd8,
  features = features_present,
  assays = DefaultAssay(seu_cd8),
  group.by = "seurat_clusters",
  slot = "data"
)

avg_mat <- avg_list[[1]]  # genes x clusters
avg_df <- data.frame(gene = rownames(avg_mat), avg_mat, check.names = FALSE)

# Order columns: gene, then cluster columns in numeric order if possible
cluster_cols <- setdiff(colnames(avg_df), "gene")
cluster_cols_sorted <- cluster_cols[order(as.numeric(cluster_cols))]
avg_df <- avg_df[, c("gene", cluster_cols_sorted)]

write.table(
  avg_df,
  file = file.path(out_dir, "CD8_cluster_avgexpr_markers.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("DONE ✅ Saved: ", file.path(out_dir, "CD8_cluster_avgexpr_markers.tsv"))

# ----------------------------
# DotPlot by CD8 clusters
# ----------------------------
p_dot <- DotPlot(
  seu_cd8,
  features = features_present,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(10, 20, 10, 10)
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "DotPlot_CD8_by_cluster.png"),
  plot = p_dot,
  width = 14, height = 6, dpi = 600
)
message("DONE ✅ Saved: ", file.path(out_dir, "DotPlot_CD8_by_cluster.png"))

# ----------------------------
# OFFICIAL CD8 naming (based on your table)
# 0 = Effector transitional
# 4 = Memory / early
# 6 = Cytotoxic effector
# ----------------------------
new_names_cd8 <- c(
  "0" = "CD8_Effector_Transitional",
  "4" = "CD8_Memory_Early",
  "6" = "CD8_Cytotoxic_Effector"
)

# Sanity: ensure all clusters are covered
missing_map <- setdiff(cd8_clusters, names(new_names_cd8))
if (length(missing_map) > 0) {
  stop("Missing name mapping for cluster(s): ", paste(missing_map, collapse = ", "),
       ". Update new_names_cd8 accordingly.")
}

# Rename idents and store as metadata column
Idents(seu_cd8) <- "seurat_clusters"
seu_cd8 <- RenameIdents(seu_cd8, new_names_cd8)
seu_cd8$CD8_Annotation <- factor(as.character(Idents(seu_cd8)),
                                 levels = unname(new_names_cd8[cd8_clusters]))

# ----------------------------
# UMAP (by official CD8_Annotation)
# Assumes UMAP already computed in object.
# If not, this will still run but DimPlot may fail.
# ----------------------------
if (!("umap" %in% names(seu_cd8@reductions))) {
  stop("UMAP not found in seu_cd8@reductions. Your T-cell object must already contain UMAP coordinates.")
}

p_umap <- DimPlot(
  seu_cd8,
  group.by = "CD8_Annotation",
  label = TRUE, repel = TRUE, label.size = 5,
  pt.size = 0.25
) +
  theme_classic(base_size = 16) +
  theme(
    axis.title = element_blank(),
    plot.title = element_blank()
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "UMAP_CD8_annotated.png"),
  plot = p_umap,
  width = 10, height = 8, dpi = 600
)
message("DONE ✅ Saved: ", file.path(out_dir, "UMAP_CD8_annotated.png"))

# ----------------------------
# Save CD8 annotated object
# ----------------------------
saveRDS(seu_cd8, file.path(out_dir, "seu_cd8_annotated.rds"))
message("DONE ✅ Saved: ", file.path(out_dir, "seu_cd8_annotated.rds"))

message("\nAll done ✅ 12_Tcell_CD8_Annotated.R completed successfully.\n")
