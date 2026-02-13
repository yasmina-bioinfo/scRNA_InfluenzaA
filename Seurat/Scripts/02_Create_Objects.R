suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(stringr)
})

# -------------------------------------------------------------------
# CONFIG
# -------------------------------------------------------------------
data_dir <- "/mnt/c/Users/yasmi/OneDrive/Desktop/ScRNA SEURAT/Influanza_Tcells/Data/geo_suppl"
out_dir  <- "Results/seurat_objects"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------
# FILES
# -------------------------------------------------------------------
tar_files <- list.files(
  data_dir,
  pattern = "EmptyDrops_CR_matrix\\.tar\\.gz$",
  full.names = TRUE
)
stopifnot(length(tar_files) > 0)

get_sample <- function(f) {
  # Example: GSM7792039_IA_1_EmptyDrops_CR_matrix.tar.gz
  m <- str_match(basename(f), "GSM\\d+_([A-Z]+_\\d+)_EmptyDrops_CR_matrix\\.tar\\.gz")
  if (is.na(m[1, 2])) stop("Nom de fichier inattendu: ", basename(f))
  m[1, 2]
}

samples <- vapply(tar_files, get_sample, character(1))

# Exclude IC
keep <- !str_detect(samples, "^IC_")
tar_files <- tar_files[keep]
samples   <- samples[keep]

message("Samples kept (IC excluded): ", paste(samples, collapse = ", "))
stopifnot(length(samples) > 0)

# -------------------------------------------------------------------
# READ CellRanger-style tar.gz (EmptyDrops output)
# Robust: finds matrix/barcodes/genes even if nested
# -------------------------------------------------------------------
read_cr_tar <- function(tar_path) {
  exdir <- tempfile("cr_")
  dir.create(exdir)
  untar(tar_path, exdir = exdir)

  # Find files recursively (some archives place them in a folder)
  mtx <- list.files(exdir, pattern = "^matrix\\.mtx$", full.names = TRUE, recursive = TRUE)
  bar <- list.files(exdir, pattern = "^barcodes\\.tsv$", full.names = TRUE, recursive = TRUE)
  gen <- list.files(exdir, pattern = "^genes\\.tsv$", full.names = TRUE, recursive = TRUE)

  if (length(mtx) == 0 || length(bar) == 0 || length(gen) == 0) {
    stop("Fichiers attendus manquants (matrix.mtx / barcodes.tsv / genes.tsv) dans: ", basename(tar_path))
  }

  mat <- readMM(mtx[1])
  # Convert to dgCMatrix (Seurat standard) to avoid class surprises
  mat <- as(mat, "dgCMatrix")

  barcodes <- readLines(bar[1])
  genes <- read.delim(gen[1], header = FALSE, stringsAsFactors = FALSE)

  gene_names <- if (ncol(genes) >= 2) genes[[2]] else genes[[1]]

  rownames(mat) <- make.unique(gene_names)
  colnames(mat) <- make.unique(barcodes)

  mat
}

# -------------------------------------------------------------------
# SAMPLE → metadata (sample / group / infection / donor_type)
# -------------------------------------------------------------------
meta_from_sample <- function(sample) {
  group <- sub("_\\d+$", "", sample) # IA, HC, PI, PHC
  infection <- ifelse(group %in% c("IA", "PI"), "infected", "healthy")
  donor_type <- ifelse(group %in% c("PHC", "PI"), "pregnant", "adult")
  list(sample = sample, group = group, infection = infection, donor_type = donor_type)
}

# -------------------------------------------------------------------
# BUILD Seurat objects
# -------------------------------------------------------------------
seu_list <- list()

for (i in seq_along(tar_files)) {
  s <- samples[i]
  message("\n>> Loading ", s)

  counts <- read_cr_tar(tar_files[i])

  obj <- CreateSeuratObject(
    counts = counts,
    project = "Influenza_Patients",
    min.cells = 3,
    min.features = 200
  )

  message("Cells: ", ncol(obj), " | Features: ", nrow(obj))
  message("Example cell name: ", colnames(obj)[1])

  # --- SAFE metadata attach: rownames MUST match cell names
  md <- meta_from_sample(s)
  cell_names <- colnames(obj)

  meta_df <- data.frame(
    sample = rep(md$sample, length(cell_names)),
    group = rep(md$group, length(cell_names)),
    infection = rep(md$infection, length(cell_names)),
    donor_type = rep(md$donor_type, length(cell_names)),
    stringsAsFactors = FALSE
  )
  rownames(meta_df) <- cell_names

  # AddMetaData will now find 100% overlap
  obj <- AddMetaData(obj, metadata = meta_df)

  # quick sanity check
  message("Meta check (sample): ", unique(obj$sample))

  seu_list[[s]] <- obj

  # optional: clean memory a bit on laptop
  rm(counts, meta_df)
  gc(verbose = FALSE)
}

# -------------------------------------------------------------------
# SAVE
# -------------------------------------------------------------------
saveRDS(seu_list, file = file.path(out_dir, "seu_list_raw.rds"))
message("\nDONE: ", length(seu_list), " Seurat objects created (IC excluded). Saved to: ", file.path(out_dir, "seu_list_raw.rds"))
