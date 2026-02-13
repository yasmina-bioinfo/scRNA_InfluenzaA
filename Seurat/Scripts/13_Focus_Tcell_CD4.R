#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

cat("=== 13_Focus_Tcell_CD4 ===\n")

# ----------------------------
# 0) Paths + checks
# ----------------------------
in_rds  <- "Results/tcells/seu_cd4_annotated.rds"
outdir  <- "Results/tcells"
out_rds <- file.path(outdir, "seu_cd4_annotated_simplified.rds")

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!file.exists(in_rds)) stop("Seurat RDS not found: ", in_rds)

# ----------------------------
# 1) Load object
# ----------------------------
seu_cd4 <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")
cat("Cells:", ncol(seu_cd4), " | Features:", nrow(seu_cd4), "\n")

# Required meta columns
required_cols <- c("seurat_clusters")
missing_cols <- setdiff(required_cols, colnames(seu_cd4@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing meta.data column(s): ", paste(missing_cols, collapse = ", "))
}

# Ensure factor
seu_cd4$seurat_clusters <- as.factor(seu_cd4$seurat_clusters)
cat("Cluster levels:", paste(levels(seu_cd4$seurat_clusters), collapse = ", "), "\n")
print(table(seu_cd4$seurat_clusters))

# Ensure UMAP exists
if (!"umap" %in% names(seu_cd4@reductions)) {
  stop("No UMAP found in object (reduction 'umap' missing). Run RunUMAP() in your reclustering step.")
}

# ----------------------------
# 2) DotPlot (CD4 panel) by cluster
# ----------------------------
features_cd4 <- c(
  # naive/memory
  "CCR7","IL7R","LEF1","TCF7","LTB",
  # cytotoxic
  "GZMB","PRF1","NKG7","GNLY",
  # Th1 axis
  "TBX21","IFNG","CXCR3",
  # Th2/Th17 (often low/absent; keep for completeness)
  "GATA3","IL4R","RORC","CCR6","IL17A",
  # Treg
  "FOXP3","IL2RA","CTLA4",
  # IFN/ISG
  "ISG15","STAT1","IFI27"
)

features_present <- features_cd4[features_cd4 %in% rownames(seu_cd4)]
features_missing <- setdiff(features_cd4, features_present)

if (length(features_present) < 6) {
  stop(
    "Too few CD4 markers found. Present: ",
    paste(features_present, collapse = ", "),
    " | Missing: ", paste(features_missing, collapse = ", ")
  )
}

p_dot <- DotPlot(
  object = seu_cd4,
  features = features_present,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(10, 20, 10, 10),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11)
  )

ggsave(
  filename = file.path(outdir, "DotPlot_CD4_by_cluster.png"),
  plot = p_dot, width = 13.5, height = 5.5, dpi = 300
)
cat("DONE ✅ Saved:", file.path(outdir, "DotPlot_CD4_by_cluster.png"), "\n")

# ----------------------------
# 3) Cluster-average table (AvgExpr)
# ----------------------------
DefaultAssay(seu_cd4) <- DefaultAssay(seu_cd4)

avg_expr <- tryCatch({
  AverageExpression(
    object = seu_cd4,
    assays = DefaultAssay(seu_cd4),
    features = features_present,
    group.by = "seurat_clusters",
    slot = "data"
  )[[DefaultAssay(seu_cd4)]]
}, error = function(e) {
  cat("Warning: AverageExpression(slot='data') failed; retrying with slot='counts'\n")
  AverageExpression(
    object = seu_cd4,
    assays = DefaultAssay(seu_cd4),
    features = features_present,
    group.by = "seurat_clusters",
    slot = "counts"
  )[[DefaultAssay(seu_cd4)]]
})

avg_df <- as.data.frame(avg_expr)
avg_df$gene <- rownames(avg_df)
avg_df <- avg_df[, c("gene", setdiff(colnames(avg_df), "gene"))]

out_avg <- file.path(outdir, "CD4_cluster_avgexpr_markers.tsv")
write.table(avg_df, out_avg, sep = "\t", quote = FALSE, row.names = FALSE)
cat("DONE ✅ Saved:", out_avg, "\n")

# ----------------------------
# 4) Simplify CD4 annotation (final states)
# ----------------------------
map_simplified <- c(
  "0" = "CD4_Naive_Memory",
  "1" = "CD4_Naive_Memory",
  "2" = "CD4_Th1_like",
  "3" = "CD4_IFN_high",
  "4" = "CD4_Th1_like",
  "5" = "CD4_IFN_high",
  "6" = "CD4_Cytotoxic_like",
  "7" = "CD4_Treg"
)

cl <- as.character(seu_cd4$seurat_clusters)
missing_keys <- setdiff(unique(cl), names(map_simplified))
if (length(missing_keys) > 0) {
  stop("Mapping missing for cluster(s): ", paste(missing_keys, collapse = ", "))
}

new_lab <- unname(map_simplified[cl])
names(new_lab) <- colnames(seu_cd4)
seu_cd4$CD4_Annotation_simplified <- factor(new_lab)

# Set biologically meaningful order (for legend & plots)
seu_cd4$CD4_Annotation_simplified <- factor(
  seu_cd4$CD4_Annotation_simplified,
  levels = c(
    "CD4_Naive_Memory",
    "CD4_Th1_like",
    "CD4_Cytotoxic_like",
    "CD4_IFN_high",
    "CD4_Treg"
  )
)

cat("CD4_Annotation_simplified levels:",
    paste(levels(seu_cd4$CD4_Annotation_simplified), collapse = ", "), "\n")

# ----------------------------
# 5) CD4 palette (fixed)
# ----------------------------
cd4_palette <- c(
  "CD4_Naive_Memory"   = "#495057",
  "CD4_Th1_like"       = "#0077b6",
  "CD4_Cytotoxic_like" = "#d00000",
  "CD4_IFN_high"       = "#ffb703",
  "CD4_Treg"           = "#9b5de5"
)

# ----------------------------
# 6) Final UMAP (annotated, no internal title)
# ----------------------------
p_umap <- DimPlot(
  object = seu_cd4,
  reduction = "umap",
  group.by = "CD4_Annotation_simplified",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  cols = cd4_palette
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12),
    plot.margin  = margin(10, 10, 10, 10)
  )

ggsave(
  file.path(outdir, "UMAP_CD4_annotated.png"),
  plot = p_umap, width = 8.5, height = 8.5, dpi = 300
)
cat("DONE ✅ Saved:", file.path(outdir, "UMAP_CD4_annotated.png"), "\n")

# ----------------------------
# 7) Save simplified object
# ----------------------------
saveRDS(seu_cd4, out_rds)
cat("DONE ✅ Saved:", out_rds, "\n")

cat("=== 13_Focus_Tcell_CD4 DONE ===\n")
