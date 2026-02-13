# This script generates pseudobulk count and metadata matrices
# for Activated-like T cells, aggregated at the patient/sample level.
# Outputs:
# - pseudobulk_counts_Activated_like_T.tsv
# - pseudobulk_meta_Activated_like_T.tsv

#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

# --- Robust paths: resolve relative to this script location ---
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", args[grep("^--file=", args)])
script_dir <- dirname(normalizePath(script_path))

project_dir <- normalizePath(file.path(script_dir, ".."))  # project root
counts_path <- file.path(project_dir, "results", "pseudobulk", "pseudobulk_counts_Activated_like_T.tsv")
meta_path   <- file.path(project_dir, "results", "pseudobulk", "pseudobulk_meta_Activated_like_T.tsv")

message("Project dir: ", project_dir)
message("Counts path: ", counts_path)
message("Meta path:   ", meta_path)

counts <- read.table(counts_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
meta   <- read.table(meta_path,   header=TRUE, row.names=1, sep="\t")

out1 <- file.path(project_dir, "results", "pseudobulk", "deseq2_Activated_like_T_IA_vs_HC.csv")
out2 <- file.path(project_dir, "results", "pseudobulk", "deseq2_Activated_like_T_PI_vs_PHC.csv")
