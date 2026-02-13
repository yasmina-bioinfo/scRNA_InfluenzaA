suppressPackageStartupMessages({
  library(Seurat)
})

# Load objects with percent.mt already computed
seu_list <- readRDS("Results/seurat_objects/seu_list_qc_added_percent_mt.rds")

dir.create("Results/qc", recursive = TRUE, showWarnings = FALSE)

filter_summary <- data.frame()

for (nm in names(seu_list)) {

  obj <- seu_list[[nm]]
  before <- ncol(obj)

  obj <- subset(
    obj,
    subset = nFeature_RNA > 200 &
             nFeature_RNA < 4000 &
             percent.mt < 10
  )

  after <- ncol(obj)

  filter_summary <- rbind(filter_summary,
    data.frame(
      sample = nm,
      cells_before = before,
      cells_after = after,
      cells_removed = before - after
    )
  )

  seu_list[[nm]] <- obj
}

# Save filtered objects
saveRDS(seu_list, "Results/seurat_objects/seu_list_filtered.rds")

# Save filtering summary
write.csv(filter_summary, "Results/qc/filtering_summary.csv", row.names = FALSE)

message("Filtering complete.")
message("Filtered objects saved to Results/seurat_objects/seu_list_filtered.rds")
message("Filtering summary saved to Results/qc/filtering_summary.csv")