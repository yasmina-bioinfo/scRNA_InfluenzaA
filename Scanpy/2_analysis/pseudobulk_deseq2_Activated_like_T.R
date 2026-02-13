#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

# --- Robust paths ---
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- dirname(normalizePath(script_path))
project_dir <- normalizePath(file.path(script_dir, ".."))

counts_path <- file.path(project_dir, "results", "pseudobulk", "pseudobulk_counts_Activated_like_T.tsv")
meta_path   <- file.path(project_dir, "results", "pseudobulk", "pseudobulk_meta_Activated_like_T.tsv")

out_adults  <- file.path(project_dir, "results", "pseudobulk", "pseudobulk_deseq2_Activated_like_T_IA_vs_HC.csv")
out_preg    <- file.path(project_dir, "results", "pseudobulk", "pseudobulk_deseq2_Activated_like_T_PI_vs_PHC.csv")

message("Project dir: ", project_dir)
message("Counts path: ", counts_path)
message("Meta path:   ", meta_path)

# --- Load inputs ---
# IMPORTANT: counts file has samples as ROWS and genes as COLUMNS (from the Python export)
counts_df <- read.table(counts_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta      <- read.table(meta_path,   header=TRUE, row.names=1, sep="\t")

# Ensure meta rownames are sample IDs
if ("sample_id" %in% colnames(meta)) rownames(meta) <- as.character(meta$sample_id)

# Keep only common samples
common <- intersect(rownames(counts_df), rownames(meta))
counts_df <- counts_df[common, , drop=FALSE]
meta <- meta[common, , drop=FALSE]

# Reorder meta to match counts_df rows
meta <- meta[rownames(counts_df), , drop=FALSE]
stopifnot(all(rownames(counts_df) == rownames(meta)))

# Transpose to DESeq2 format: genes (rows) x samples (cols)
counts_mat <- t(as.matrix(counts_df))
mode(counts_mat) <- "integer"

# Quick sanity check: should NOT be all zeros
message("Nonzero entries in counts_mat: ", sum(counts_mat > 0))
message("Library sizes (samples):")
print(colSums(counts_mat))

run_deseq <- function(conds, ref_level, out_csv) {
  keep <- meta$condition %in% conds
  meta_sub <- meta[keep, , drop=FALSE]

  # subset counts columns by sample IDs
  counts_sub <- counts_mat[, rownames(meta_sub), drop=FALSE]

  meta_sub$condition <- factor(meta_sub$condition, levels=c(ref_level, setdiff(conds, ref_level)))

  dds <- DESeqDataSetFromMatrix(
    countData = counts_sub,
    colData = meta_sub,
    design = ~ condition
  )

  dds <- DESeq(dds)
  res <- results(dds)
  res <- res[order(res$padj), ]
  write.csv(as.data.frame(res), out_csv)
  message("Saved: ", out_csv)
}

# Adults
run_deseq(c("HC","IA"), "HC", out_adults)

# Pregnancy
run_deseq(c("PHC","PI"), "PHC", out_preg)
