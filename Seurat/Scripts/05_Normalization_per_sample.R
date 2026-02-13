suppressPackageStartupMessages({
  library(Seurat)
})

seu_list <- readRDS("Results/seurat_objects/seu_list_filtered.rds")
dir.create("Results/seurat_objects", recursive = TRUE, showWarnings = FALSE)

seu_list <- readRDS(in_rds)

nfeatures_hvg <- 2000
npcs <- 30  # Reduced number for Laptop

for (nm in names(seu_list)) {
  message("\n>> Normalisation + HVG (léger) pour: ", nm)
  obj <- seu_list[[nm]]
  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = nfeatures_hvg,
    verbose = FALSE
  )

  # Scale only HVG for PCA
  obj <- ScaleData(
    obj,
    features = VariableFeatures(obj),
    verbose = FALSE
  )

  # PCA on HVG
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = npcs,
    verbose = FALSE
  )

  seu_list[[nm]] <- obj

  # Memory management
  gc(verbose = FALSE)
}

saveRDS(seu_list, out_rds)
message("\nDONE: normalization + HVG + PCA saved -> ", out_rds)
