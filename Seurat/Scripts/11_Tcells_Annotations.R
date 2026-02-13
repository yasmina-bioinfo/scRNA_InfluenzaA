# scripts/11_Tcells_Annotation.R
# Purpose: Load reclustered T-cell subset, add named T-cell annotations,
#          generate an annotated UMAP, and save outputs.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# 0) Input / Output
# -----------------------------
in_rds  <- "Results/tcells/seu_tcells_subset_reclustered.rds"
out_dir <- "Results/tcells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(in_rds)) stop("Input RDS not found: ", in_rds)

# -----------------------------
# 1) Load object
# -----------------------------
seu_t <- readRDS(in_rds)

# -----------------------------
# 2) Checks
# -----------------------------
if (!"seurat_clusters" %in% colnames(seu_t@meta.data)) {
  stop("Column 'seurat_clusters' not found in meta.data.")
}
if (!"umap" %in% names(seu_t@reductions)) {
  stop("UMAP reduction not found in object.")
}

seu_t$seurat_clusters <- as.factor(seu_t$seurat_clusters)

# -----------------------------
# 3) Cluster -> named annotation
# -----------------------------
tcell_labels <- c(
  "0"  = "CD8_Effector",
  "1"  = "CD4_Naive_TCM",
  "2"  = "CD4_Effector_Memory",
  "3"  = "ISG_high_T",
  "4"  = "CD8_Memory",
  "5"  = "CD4_Activated_IFN",
  "6"  = "CD8_Terminal_Cytotoxic",
  "7"  = "GammaDelta_T",
  "8"  = "Atypical_T",
  "9"  = "Myeloid_like?",
  "10" = "MAIT_like",
  "11" = "Myeloid_like?"
)

seu_t$Tcell_Annotation <- unname(tcell_labels[as.character(seu_t$seurat_clusters)])

# Validate: no missing labels
if (any(is.na(seu_t$Tcell_Annotation))) {
  missing <- sort(unique(as.character(seu_t$seurat_clusters[is.na(seu_t$Tcell_Annotation)])))
  stop("Missing Tcell_Annotation for clusters: ", paste(missing, collapse = ", "),
       "\nUpdate tcell_labels mapping.")
}

# Save a summary table
summary_tbl <- as.data.frame(table(seu_t$seurat_clusters, seu_t$Tcell_Annotation))
colnames(summary_tbl) <- c("seurat_clusters", "Tcell_Annotation", "n_cells")

write.table(summary_tbl,
            file = file.path(out_dir, "tcell_annotation_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\n=== Tcell_Annotation summary ===\n")
print(summary_tbl)

# -----------------------------
# 4) UMAP plots
# -----------------------------
p_annot <- DimPlot(
  seu_t,
  reduction = "umap",
  group.by = "Tcell_Annotation",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.2
) +
  theme_classic(base_size = 13) +
  theme(
    legend.title = element_blank(),
    legend.text  = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    axis.title = element_text(size = 12),
    axis.text  = element_text(size = 10),
    plot.title = element_blank()
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4)
    )
  )

ggplot2::ggsave(
  filename = file.path(out_dir, "UMAP_Tcells_annotated.png"),
  plot = p_annot,
  width = 7,
  height = 7,
  units = "in",
  dpi = 600
)

# -----------------------------
# 5) Save updated object
# -----------------------------
saveRDS(seu_t, file = file.path(out_dir, "seu_tcells_annotated.rds"))

message("\nDONE ✅ Saved:\n",
        "- Results/tcells/tcell_annotation_summary.tsv\n",
        "- Results/tcells/UMAP_Tcells_annotated.png\n",
        "- Results/tcells/seu_tcells_annotated.rds\n")

# -----------------------------
# 6) DotPlot global (T cells)
# -----------------------------
DefaultAssay(seu_t) <- "RNA"

# Panel global (T cells): lineage + CD4/CD8 + cytotoxic + naive/memory + ISG + activation/exhaustion + special subsets
features_global <- c(
  # T cell core
  "TRAC", "CD3D", "CD3E",
  # CD4 / naive-memory
  "CD4", "IL7R", "CCR7", "LEF1", "LTB",
  # CD8 / cytotoxic
  "CD8A", "CD8B", "NKG7", "GNLY", "GZMB", "PRF1", "FGFBP2",
  # Exhaustion / activation
  "LAG3", "PDCD1", "TIGIT", "TOX", "NR4A1", "NR4A2",
  # ISG / interferon
  "ISG15", "IFIT2", "IFI27", "OAS1", "OASL",
  # Gamma-delta
  "TRDC", "TRGC1", "TRGC2",
  # MAIT-like
  "KLRB1", "TRAV1-2", "SLC4A10"
)

# Keep only genes present in dataset (avoid errors)
features_global_present <- features_global[features_global %in% rownames(seu_t)]
if (length(features_global_present) < 5) stop("Too few markers found in object. Check gene symbols.")

p <- DotPlot(
  seu_t,
  features = features_present,
  group.by = "Tcell_Annotation",
  dot.scale = 6
) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 11),
    plot.margin = margin(10, 20, 10, 10)  # un peu d'air
  )

ggsave(
  "Results/tcells/DotPlot_Tcells_Global.png",
  plot = p,
  width = 16,
  height = 6,
  units = "in",
  dpi = 600
)

