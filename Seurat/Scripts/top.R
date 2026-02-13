
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})
suppressPackageStartupMessages({
  library(dplyr)
})

# Recharge les markers déjà sauvés
markers <- read.delim("Results/markers/markers_all_clusters_downsample.tsv",
                      stringsAsFactors = FALSE, check.names = FALSE)

# ---- Condensed summary (one line per cluster) ----
summary_clusters <- markers %>%
  mutate(gene = if ("gene" %in% colnames(markers)) gene else rownames(markers)) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  summarise(top_genes = paste(gene, collapse = ", "), .groups = "drop") %>%
  arrange(as.numeric(as.character(cluster)))

write.table(summary_clusters,
            "Results/annotation/top20_markers_per_cluster.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

message("\n=== SUMMARY (Top 20 genes per cluster) ===\n")
print(summary_clusters, n = Inf)
