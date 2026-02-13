# scRNA-seq Influenza Patients — T Cell–Focused Analysis

## Project Summary
This project presents an end-to-end single-cell RNA-seq (scRNA-seq) analysis of peripheral blood immune cells from influenza patients, with a focused investigation of **T cell transcriptional responses** to infection.

The analysis combines:
- cell-level scRNA-seq exploration,
- targeted T cell subpopulation analysis,
- patient-level **pseudobulk aggregation**, and
- **differential expression analysis using DESeq2**.

The workflow is designed to reflect **realistic research conditions**, including computational constraints, data quality considerations, and biologically driven analytical choices.

---

## Biological Question
How does influenza infection alter **T cell transcriptional states** at the patient level, and are these changes consistent across:
- infected vs healthy adults,
- infected vs healthy pregnant individuals?

---

## Dataset
- **Source:** GEO accession **GSE243629**
- **Samples:** Peripheral blood mononuclear cells (PBMCs)
- **Conditions:**
  - HC: Healthy adults
  - IA: Influenza-infected adults
  - PHC: Healthy pregnant controls
  - PI: Influenza-infected pregnant patients

Raw sequencing data are publicly available and are not stored in this repository.

---

## Analysis Overview

### 1. Quality Control and Preprocessing
- Cell- and gene-level quality control.
- Removal of low-quality cells.
- Normalization and selection of highly variable genes.
- Dimensionality reduction, neighborhood graph construction, and UMAP embedding.

### 2. PBMC Clustering and Annotation
- Unsupervised clustering of PBMCs.
- Marker-based annotation of major immune cell types.
- Quantification of cell-type proportions across conditions.

### 3. T Cell–Focused Analysis
- Extraction of T cells from the annotated PBMC dataset
- Subclustering and annotation of major T cell functional states:
  - CD4 naive/memory-like
  - Cytotoxic T cells
  - Activated-like T cells
  - Cycling T cells
- Visualization using UMAP and summary plots.

### 4. Pseudobulk Strategy
- Activated-like T cells were selected for downstream analysis to ensure sufficient counts per patient.
- Raw gene counts were aggregated per patient (pseudobulk).
- Separate comparisons were performed:
  - **Adults:** IA vs HC
  - **Pregnancy:** PI vs PHC

### 5. Differential Expression Analysis
- Differential expression analysis performed using **DESeq2**
- Results exported as tables for downstream interpretation.
- This approach enables robust patient-level inference while preserving upstream single-cell resolution.

---

## Data availability

Large intermediate objects (e.g. `.h5ad` files, raw matrices) are not included in this repository due to size constraints.
The repository focuses on reproducible analysis code, figures, and interpretation.

---

## Key Results
- Influenza infection is associated with an increased proportion of Activated-like T cells across patient groups.
- Differential expression analysis of pseudobulked Activated-like T cells highlights activation-associated transcriptional programs in influenza-infected individuals.
- Distinct transcriptional patterns are observed between adult and pregnant cohorts, suggesting condition-specific modulation of T cell activation.

The `results/` directory contains:
- `figures/`: final figures used in the report
- other subfolders: placeholders documenting the analysis structure (not versioned due to file size).

A dedicated `Figures_only.ipynb` notebook is provided to regenerate publication-quality figures from final AnnData objects, without re-running the full pipeline.

Overall, the analysis suggests that influenza reshapes transcriptional programs within existing T cell populations rather than creating new cell types.
---

## Repository Structure

    scRNA_Influenza_Patients/
    ├── 1_qc/                     # Quality control notebooks
    ├── 2_analysis/               # Main analysis notebooks and scripts
    │   ├── pseudobulk_build_Activated_like_T.py
    │   └── pseudobulk_deseq2_Activated_like_T.R
    ├── data_raw/                 # Raw data (not tracked)
    ├── results/
    │   ├── figures/              # UMAPs and summary plots
    │   ├── pseudobulk/           # Pseudobulk matrices and DESeq2 results
    │   └── README.md
    ├── .gitignore
    └── README.md                 # This file

Large intermediate files (e.g. `.h5ad`, raw sequencing archives) are excluded from version control.

---

## Notes and Limitations
- Pseudobulk analysis was restricted to Activated-like T cells to ensure sufficient patient-level counts.
- CD4/CD8 separation was not applied at the pseudobulk stage by design.
- Parameter choices reflect computational limitations typical of personal computing environments.
- All analytical decisions are explicitly documented in the notebooks.

---

## Purpose of This Project
This project provides a reproducible scRNA-seq case study and a biologically grounded analysis of T cell responses to influenza infection.

---

## Citation

Zhang Y. et al. (2023). *A single-cell atlas of the peripheral immune response in patients with influenza A virus infection*.  
GEO accession: GSE243629.

