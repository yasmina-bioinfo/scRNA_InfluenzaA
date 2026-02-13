# ===============================
# Temporary script – DotPlot T cells
# ===============================

library(Seurat)
library(ggplot2)

# 1) Load T cell object
seu_t <- readRDS("Results/tcells/seu_tcells_annotated.rds")

# 2) Marker panel
features_global <- c(
  "TRAC", "CD3D", "CD3E",
  "CD4", "IL7R", "CCR7", "LEF1", "LTB",
  "CD8A", "CD8B", "NKG7", "GNLY", "GZMB", "PRF1", "FGFBP2",
  "LAG3", "PDCD1", "TIGIT", "TOX",
  "ISG15", "IFIT2", "IFI27", "OAS1", "OASL",
  "TRDC", "TRGC1", "TRGC2",
  "KLRB1", "TRAV1-2", "SLC4A10"
)

features_present <- features_global[features_global %in% rownames(seu_t)]

if (length(features_present) < 5) {
  stop("Too few markers found in object. Check gene symbols.")
}

# 3) DotPlot
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

# 4) Save
dir.create("Results/tcells", showWarnings = FALSE)

ggsave(
  "Results/tcells/DotPlot_Tcells_Global.png",
  plot = p,
  width = 16,
  height = 6,
  units = "in",
  dpi = 600
)

