suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

# -----------------------------
# Paths (project root)
# -----------------------------
in_rds  <- "Results/seurat_objects/seu_list_raw.rds"
out_rds <- "Results/seurat_objects/seu_list_qc_added_percent_mt.rds"

dir.create("Results/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/figures/qc", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/seurat_objects", recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Load
# -----------------------------
seu_list <- readRDS(in_rds)

qc_summary <- data.frame()

# Helper: stable violin plot from meta.data (no Seurat plotting)
plot_violin_qc <- function(df, sample_name) {
  df$metric <- factor(df$metric, levels = c("nFeature_RNA", "nCount_RNA", "percent.mt"))

  ggplot(df, aes(x = metric, y = value)) +
    geom_violin(trim = TRUE) +
    geom_jitter(width = 0.15, size = 0.15, alpha = 0.25) +
    theme_bw() +
    labs(title = paste0("QC violins - ", sample_name), x = NULL, y = NULL)
}

plot_scatter <- function(df, x, y, sample_name, title_suffix) {
  ggplot(df, aes_string(x = x, y = y)) +
    geom_point(size = 0.3, alpha = 0.35) +
    theme_bw() +
    labs(title = paste0(title_suffix, " - ", sample_name), x = x, y = y)
}

# -----------------------------
# QC per sample
# -----------------------------
for (nm in names(seu_list)) {
  obj <- seu_list[[nm]]

  # Compute percent.mt (human mito genes)
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")

  md <- obj@meta.data

  # Summary stats
  qc_summary <- rbind(qc_summary, data.frame(
    sample = nm,
    cells = ncol(obj),
    median_nFeature = median(md$nFeature_RNA),
    median_nCount = median(md$nCount_RNA),
    median_percent_mt = median(md$percent.mt),
    stringsAsFactors = FALSE
  ))

  # Prepare long format for violins
  df_long <- rbind(
    data.frame(metric = "nFeature_RNA", value = md$nFeature_RNA),
    data.frame(metric = "nCount_RNA",   value = md$nCount_RNA),
    data.frame(metric = "percent.mt",   value = md$percent.mt)
  )

  # Violin QC
  p_vln <- plot_violin_qc(df_long, nm)
  ggsave(
    filename = file.path("Results/figures/qc", paste0("QC_violin_", nm, ".pdf")),
    plot = p_vln, width = 10, height = 4
  )

  # Scatter QC
  p_sc1 <- plot_scatter(md, "nCount_RNA", "percent.mt", nm, "nCount vs percent.mt")
  ggsave(
    filename = file.path("Results/figures/qc", paste0("QC_scatter_nCount_vs_mt_", nm, ".pdf")),
    plot = p_sc1, width = 5.5, height = 4.5
  )

  p_sc2 <- plot_scatter(md, "nCount_RNA", "nFeature_RNA", nm, "nCount vs nFeature")
  ggsave(
    filename = file.path("Results/figures/qc", paste0("QC_scatter_nCount_vs_nFeature_", nm, ".pdf")),
    plot = p_sc2, width = 5.5, height = 4.5
  )

  seu_list[[nm]] <- obj
}

# Save outputs
write.csv(qc_summary, "Results/qc/qc_summary_by_sample.csv", row.names = FALSE)
saveRDS(seu_list, out_rds)

message("DONE: QC metrics computed + figures saved in Results/figures/qc/")
message("QC summary: Results/qc/qc_summary_by_sample.csv")
message("Updated objects: ", out_rds)
