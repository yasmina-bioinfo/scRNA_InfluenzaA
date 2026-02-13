# -----------------------------
# DotPlot CD8 (fine, by seurat_clusters)
# -----------------------------
library(Seurat)
library(ggplot2)

seu_t <- readRDS("Results/tcells/seu_tcells_annotated.rds")
DefaultAssay(seu_t) <- "RNA"

features_cd8 <- c(
  "CD8A","CD8B",
  "NKG7","GZMB","PRF1","GNLY",
  "PDCD1","LAG3","TIGIT","TOX",
  "CCR7","IL7R",
  "IFNG","STAT1","ISG15",
  "MKI67"
)

features_present <- features_cd8[features_cd8 %in% rownames(seu_t)]
if (length(features_present) < 10) stop("Too few CD8 markers found. Check gene symbols.")

p_cd8 <- DotPlot(
  seu_t,
  features = features_present,
  group.by = "seurat_clusters",
  dot.scale = 6
) +
  theme_classic(base_size = 12) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    plot.margin = margin(10, 20, 10, 10)
  )

ggsave(
  "Results/tcells/DotPlot_CD8_by_cluster.png",
  plot = p_cd8,
  width = 14,
  height = 6,
  units = "in",
  dpi = 600
)