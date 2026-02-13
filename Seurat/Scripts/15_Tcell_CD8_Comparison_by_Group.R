#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

cat("=== 16_Tcell_CD8_Comparison_by_Group ===\n")

# ----------------------------
# 0) Paths
# ----------------------------
in_rds <- "Results/tcells/seu_cd8_annotated.rds"   # adapte si besoin
outdir <- "Results/tcells"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
if (!file.exists(in_rds)) stop("RDS not found: ", in_rds)

# ----------------------------
# 1) Load + checks
# ----------------------------
seu_cd8 <- readRDS(in_rds)
cat("Loaded:", in_rds, "\n")
cat("Cells:", ncol(seu_cd8), " | Features:", nrow(seu_cd8), "\n")

req_cols <- c("group", "CD8_Annotation")
missing <- setdiff(req_cols, colnames(seu_cd8@meta.data))
if (length(missing) > 0) stop("Missing meta columns: ", paste(missing, collapse = ", "))
if (!"umap" %in% names(seu_cd8@reductions)) stop("UMAP not found in reductions.")

# Build ordered 4-group condition column
desired_levels <- c("HC", "IA", "PHC", "PI")
seu_cd8$cond4 <- as.character(seu_cd8$group)

missing_levels <- setdiff(desired_levels, unique(seu_cd8$cond4))
if (length(missing_levels) > 0) {
  stop("Missing expected group levels: ", paste(missing_levels, collapse = ", "),
       " | Present: ", paste(sort(unique(seu_cd8$cond4)), collapse = ", "))
}
seu_cd8$cond4 <- factor(seu_cd8$cond4, levels = desired_levels)

cat("cond4 levels:", paste(levels(seu_cd8$cond4), collapse = ", "), "\n")
print(table(seu_cd8$cond4))

# Enforce CD8 state order (legend/barplot)
seu_cd8$CD8_Annotation <- factor(
  seu_cd8$CD8_Annotation,
  levels = c("CD8_Memory_Early", "CD8_Effector_Transitional", "CD8_Cytotoxic_Effector")
)

# CD8 palette (your style: blue family for T/NK)
cd8_palette <- c(
  "CD8_Memory_Early"          = "#a6cee3",  # bleu clair
  "CD8_Effector_Transitional" = "#1f78b4",  # bleu moyen
  "CD8_Cytotoxic_Effector"    = "#03045e"   # bleu foncé
)

# ----------------------------
# 2) Barplot (CD8 states by group)
# ----------------------------
df_cd8 <- seu_cd8@meta.data %>%
  mutate(cond4 = factor(as.character(group), levels = desired_levels)) %>%
  group_by(cond4, CD8_Annotation) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cond4) %>%
  mutate(freq = n / sum(n))

p_bar <- ggplot(df_cd8, aes(x = cond4, y = freq, fill = CD8_Annotation)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12)
  ) +
  labs(x = NULL, y = "Proportion") +
  scale_fill_manual(values = cd8_palette)

ggsave(file.path(outdir, "Barplot_CD8_states_HC_IA_PHC_PI.png"),
       p_bar, width = 7.5, height = 6, dpi = 300)
cat("DONE ✅ Saved: Barplot_CD8_states_HC_IA_PHC_PI.png\n")

# ----------------------------
# 3) UMAP (by group) + UMAP (CD8 states)
# ----------------------------
p_umap_group <- DimPlot(
  seu_cd8, reduction = "umap", group.by = "cond4"
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12)
  )

ggsave(file.path(outdir, "UMAP_CD8_HC_IA_PHC_PI.png"),
       p_umap_group, width = 8.5, height = 8.5, dpi = 300)
cat("DONE ✅ Saved: UMAP_CD8_HC_IA_PHC_PI.png\n")

p_umap_states <- DimPlot(
  seu_cd8,
  reduction = "umap",
  group.by = "CD8_Annotation",
  label = TRUE,
  repel = TRUE,
  label.size = 5,
  cols = cd8_palette
) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_blank(),
    axis.title = element_text(size = 16),
    axis.text  = element_text(size = 12),
    legend.title = element_blank(),
    legend.text  = element_text(size = 12)
  )

ggsave(file.path(outdir, "UMAP_CD8_annotated.png"),
       p_umap_states, width = 8.5, height = 8.5, dpi = 300)
cat("DONE ✅ Saved: UMAP_CD8_annotated.png\n")

# ----------------------------
# 4) DotPlot (markers by group)
# ----------------------------
features_dot <- c(
  "CCR7","IL7R","TCF7",                # memory-like
  "GZMB","PRF1","NKG7","GNLY",         # cytotoxic
  "ISG15","STAT1",                     # IFN
  "PDCD1","LAG3","TOX"                 # exhaustion (often low in acute infection)
)
features_present <- features_dot[features_dot %in% rownames(seu_cd8)]
if (length(features_present) < 6) stop("Too few marker genes found for DotPlot.")

p_dot <- DotPlot(
  seu_cd8,
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

ggsave(file.path(outdir, "DotPlot_CD8_HC_IA_PHC_PI.png"),
       p_dot, width = 12.5, height = 4.8, dpi = 300)
cat("DONE ✅ Saved: DotPlot_CD8_HC_IA_PHC_PI.png\n")

# ----------------------------
# 5) FeaturePlot (key genes)
# ----------------------------
features_fp <- c("IL7R","CCR7","TCF7","GZMB","PRF1","NKG7","ISG15","STAT1","PDCD1","TOX")
features_fp <- features_fp[features_fp %in% rownames(seu_cd8)]
if (length(features_fp) == 0) stop("No FeaturePlot genes found in object.")

# A) Global
p_fp_global <- FeaturePlot(
  seu_cd8,
  features = features_fp,
  reduction = "umap",
  ncol = 5
)

ggsave(file.path(outdir, "FeaturePlot_CD8_global.png"),
       p_fp_global, width = 16, height = 8.5, dpi = 300)
cat("DONE ✅ Saved: FeaturePlot_CD8_global.png\n")

# B) Split by group (keep it readable)
features_fp_split <- intersect(c("IL7R","GZMB","ISG15","PDCD1"), rownames(seu_cd8))
if (length(features_fp_split) > 0) {
  p_fp_split <- FeaturePlot(
    seu_cd8,
    features = features_fp_split,
    reduction = "umap",
    split.by = "cond4",
    ncol = 4
  )
  ggsave(file.path(outdir, "FeaturePlot_CD8_split_by_group.png"),
         p_fp_split, width = 18, height = 10, dpi = 300)
  cat("DONE ✅ Saved: FeaturePlot_CD8_split_by_group.png\n")
} else {
  cat("Note: split FeaturePlot skipped (missing IL7R/GZMB/ISG15/PDCD1).\n")
}

cat("=== 16_Tcell_CD8_Comparison_by_Group DONE ===\n")
