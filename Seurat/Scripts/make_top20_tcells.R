suppressPackageStartupMessages({
  library(dplyr)
})

m <- read.delim("Results/tcells/markers_tcells_all_clusters.tsv", sep = "\t")

cluster_col <- if ("cluster" %in% colnames(m)) {
  "cluster"
} else if ("ident" %in% colnames(m)) {
  "ident"
} else {
  stop("No cluster/ident column in markers file.")
}

gene_col <- if ("gene" %in% colnames(m)) {
  "gene"
} else {
  stop("No gene column in markers file.")
}

top <- m %>%
  group_by(.data[[cluster_col]]) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  summarise(top_genes = paste(.data[[gene_col]], collapse = ", "),
            .groups = "drop") %>%
  rename(cluster = .data[[cluster_col]]) %>%
  arrange(as.numeric(as.character(cluster)))

write.table(
  top,
  "Results/tcells/top20_markers_per_tcell_cluster.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Top20 file created successfully.\n")