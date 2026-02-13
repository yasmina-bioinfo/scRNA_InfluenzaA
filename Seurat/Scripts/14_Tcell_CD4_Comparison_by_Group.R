#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

cat("=== 15_Tcell_CD4_Comparison_by_Group ===\n")

# ----------------------------
# 0) Paths
# ----------------------------
in_rds <- "Results/tcells/seu_cd4_annotated_simplified.rds"
outdir <- "Results/tcells"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!file.exists(in_rds)) stop("RDS not found: ", in_rds)

# ----------------------------
# 1) Load + checks
# ----------------------------
seu_cd4 <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")
cat("Cells:", ncol(seu_cd4), " | Features:", nrow(seu_cd4), "\n")

req_cols <- c("group", "CD4_Annotation_simplified")
missing <- setdiff(req_cols, colnames(seu_cd4@meta.data))
if (length(missing) > 0) stop("Missing meta columns: ", paste(missing, collapse = ", "))
if (!"umap" %in% names(seu_cd4@reductions)) stop("UMAP not found in reductions.")

# Build ordered 4-group condition column
desired_levels <- c("HC", "IA", "PHC", "PI")
seu_cd4$cond4 <- as.character(seu_cd4$group)

missing_levels <- setdiff(desired_levels, unique(seu_cd4$cond4))
if (length(missing_levels) > 0) {
  stop("Missing expected group levels: ", paste(missing_levels, collapse = ", "),
       " | Present: ", paste(sort(unique(seu_cd4$cond4)), collapse = ", "))
}
seu_cd4$cond4 <- factor(seu_cd4$cond4, levels = desired_levels)

cat("cond4 levels:", paste(levels(seu_cd4$cond4), collapse = ", "), "\n")
print(table(seu_cd4$cond4))

# (Optional but recommended) enforce CD4 state order (legend/barplot)
seu_cd4$CD4_Annotation_simplified <- factor(
  seu_cd4$CD4_Annotation_simplified,
  levels = c("CD4_Naive_Memory","CD4_Th1_like","CD4_Cytotoxic_like","CD4_IFN_high","CD4_Treg")
)

# CD4 palette (your fixed palette)
cd4_palette <- c(
  "CD4_Naive_Memory"   = "#495057",
  "CD4_Th1_like"       = "#0077b6",
  "CD4_Cytotoxic_like" = "#d00000",
  "CD4_IFN_high"       = "#ffb703",
  "CD4_Treg"           = "#9b5de5"
)

# ----------------------------
# 2) Barplot (CD4 states by group)
# ----------------------------
df_cd4 <- seu_cd4@meta.data %>%
  mutate(cond4 = factor(as.character(group), levels = desired_levels)) %>%
  group_by(cond4, CD4_Annotation_simplified) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cond4) %>%
  mutate(freq = n / sum(n))

p_bar <- ggplot(df_cd4, aes(x = cond4, y = freq, fill = CD4_Annotation_simplified)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12)
  ) +
  labs(x = NULL, y = "Proportion") +
  scale_fill_manual(values = cd4_palette)

ggsave(file.path(outdir, "Barplot_CD4_states_HC_IA_PHC_PI.png"),
       p_bar, width = 7.5, height = 6, dpi = 300)
cat("DONE ✅ Saved: Barplot_CD4_states_HC_IA_PHC_PI.png\n")

# ----------------------------
# 3) UMAP (by group) + UMAP (CD4 states)
# ----------------------------
p_umap_group <- DimPlot(
  seu_cd4, reduction = "umap", group.by = "cond4"
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12)
  )

ggsave(file.path(outdir, "UMAP_CD4_by_group_HC_IA_PHC_PI.png"),
       p_umap_group, width = 8.5, height = 8.5, dpi = 300)
cat("DONE ✅ Saved: UMAP_CD4_by_group_HC_IA_PHC_PI.png\n")

p_umap_states <- DimPlot(
  seu_cd4,
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
    legend.text  = element_text(size = 12)
  )

ggsave(file.path(outdir, "UMAP_CD4_annotated.png"),
       p_umap_states, width = 8.5, height = 8.5, dpi = 300)
cat("DONE ✅ Saved: UMAP_CD4_annotated.png\n")

# ----------------------------
# 4) DotPlot (markers by group)
# ----------------------------
features_dot <- c(
  "CCR7","IL7R","TCF7",          # naive/memory
  "CXCR3","IFNG",                # Th1
  "GZMB","PRF1","NKG7","GNLY",   # cytotoxic
  "ISG15","STAT1",               # IFN
  "FOXP3","CTLA4","IL2RA"        # Treg
)
features_present <- features_dot[features_dot %in% rownames(seu_cd4)]
if (length(features_present) < 6) stop("Too few marker genes found for DotPlot.")

p_dot <- DotPlot(
  seu_cd4,
  features = features_present,
  group.by = "cond4",
  dot.scale = 6
) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 11)
  )

ggsave(file.path(outdir, "DotPlot_CD4_by_group_HC_IA_PHC_PI.png"),
       p_dot, width = 11.5, height = 4.8, dpi = 300)
cat("DONE ✅ Saved: DotPlot_CD4_by_group_HC_IA_PHC_PI.png\n")

# ----------------------------
# 5) FeaturePlot (key genes, split by group)
# ----------------------------
features_fp <- c("IL7R","CCR7","CXCR3","IFNG","ISG15","STAT1","GZMB","NKG7","FOXP3","CTLA4")
features_fp <- features_fp[features_fp %in% rownames(seu_cd4)]
if (length(features_fp) == 0) stop("No FeaturePlot genes found in object.")

# A) global FeaturePlot (all cells)
p_fp_global <- FeaturePlot(
  seu_cd4,
  features = features_fp,
  reduction = "umap",
  ncol = 5
)

ggsave(file.path(outdir, "FeaturePlot_CD4_keygenes_global.png"),
       p_fp_global, width = 16, height = 8.5, dpi = 300)
cat("DONE ✅ Saved: FeaturePlot_CD4_keygenes_global.png\n")

# B) split by group (may be large; keep only 4 genes for readability)
features_fp_split <- intersect(c("IL7R","CXCR3","ISG15","FOXP3"), rownames(seu_cd4))
if (length(features_fp_split) > 0) {
  p_fp_split <- FeaturePlot(
    seu_cd4,
    features = features_fp_split,
    reduction = "umap",
    split.by = "cond4",
    ncol = 4
  )
  ggsave(file.path(outdir, "FeaturePlot_CD4_split_by_group.png"),
         p_fp_split, width = 18, height = 10, dpi = 300)
  cat("DONE ✅ Saved: FeaturePlot_CD4_split_by_group.png\n")
} else {
  cat("Note: split FeaturePlot skipped (missing IL7R/CXCR3/ISG15/FOXP3).\n")
}

cat("=== 15_Tcell_CD4_Comparison_by_Group DONE ===\n")
