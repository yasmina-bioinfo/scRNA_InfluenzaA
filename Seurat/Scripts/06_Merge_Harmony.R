suppressPackageStartupMessages({
  library(Seurat)
})

# Harmony package
# install.packages("harmony")  # Do it only if not already installed
suppressPackageStartupMessages(library(harmony))

dir.create("Results/seurat_objects", recursive = TRUE, showWarnings = FALSE)

seu_list <- readRDS("Results/seurat_objects/seu_list_filtered.rds")

# make sure cell names are unique across samples by adding sample ID as prefix
for (nm in names(seu_list)) {
  seu_list[[nm]] <- RenameCells(seu_list[[nm]], add.cell.id = nm)
  DefaultAssay(seu_list[[nm]]) <- "RNA"
}

# Merge (let's use the first sample as reference and merge the rest to it)
seu <- merge(seu_list[[1]], y = seu_list[-1], project = "Influenza_Patients")
rm(seu_list); gc(FALSE)

# Standard workflow
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 30, verbose = FALSE)

# Batch correction with Harmony using sample ID
seu <- RunHarmony(seu, group.by.vars = "sample", dims.use = 1:30, verbose = FALSE)

# UMAP/clustering on Harmony embeddings
seu <- RunUMAP(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- FindNeighbors(seu, reduction = "harmony", dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.5, verbose = FALSE)

saveRDS(seu, "Results/seurat_objects/seu_harmony_merged.rds")
message("DONE: Harmony-corrected merged object saved -> Results/seurat_objects/seu_harmony_merged.rds")
