## Influenza A Infection — Single-Cell Immune Landscape and T Cell Functional Remodeling
---
 
### Project Overview

This project presents an end-to-end single-cell RNA-seq (scRNA-seq) analysis of peripheral blood immune cells from influenza A–infected patients, with a focused investigation of T cell functional remodeling during acute infection.

The analysis was independently implemented in both Python (Scanpy) and R (Seurat) to ensure methodological robustness and cross-framework reproducibility.

The objective was to move from:

- global immune landscape characterization 
to
-refined T cell functional state analysis 
to
- condition-dependent remodeling across adult and pregnant cohorts.

---

### Biological Question

How does acute influenza infection reshape immune cell composition and T cell transcriptional states across:
- IA vs HC (infected vs healthy adults)
- PI vs PHC (infected vs healthy pregnant individuals)?

Does infection induce new T cell populations, or remodel existing functional programs?

--- 

### Dataset

- GEO accession: GSE243629
- PBMC single-cell RNA-seq
- Conditions:
    - HC — Healthy adults
    - IA — Influenza-infected adults
    - PHC — Healthy pregnant controls
    - PI — Influenza-infected pregnant patients

Raw sequencing files are not stored in this repository due to size constraints.

--- 

### Analysis Strategy

**1. Global PBMC Landscape**

- Quality control and filtering
- Normalization and dimensionality reduction
- Batch correction (Harmony in Seurat pipeline)
- Unsupervised clustering
- Broad lineage annotation (T/NK, B/Plasma, Myeloid, Neutrophils, Others)

This step established the immune context before focusing on T cells.

**2. T Cell Subsetting and Functional Annotation**

T cells were extracted and subclustered to characterize functional states.

Major identified states include:
- Naive / central memory (IL7R, CCR7, TCF7)
- Cytotoxic effector (GZMB, PRF1, NKG7, GNLY)
- IFN-responsive subsets (ISG15, STAT1)
- Th1-like CD4 populations (CXCR3, IFNG)
- Regulatory T cells (FOXP3, CTLA4)
- γδ T cells (TRDC, TRGC1/2)
- MAIT-like populations (KLRB1 + TRAV1-2)

Across both CD4 and CD8 compartments, infection was associated with:
- Expansion of cytotoxic programs
- Increased interferon-responsive states
- No dominant terminal exhaustion signature

These observations are consistent with an acute antiviral immune response rather than chronic dysfunction.

--- 

**4. Pseudobulk Patient-Level Analysis (Scanpy Pipeline)**

To enable patient-level inference:
- Activated-like T cells were pseudobulk aggregated per individual
- Separate comparisons:
    - IA vs HC
    - PI vs PHC
- Differential expression performed using DESeq2

This approach complements single-cell resolution with statistically robust patient-level comparisons.

---

### Cross-Scale Interpretation

Together, the analysis suggests that acute influenza infection:
- Redistributes immune cell states rather than creating new lineages
- Promotes coordinated CD4 helper polarization and CD8 cytotoxic differentiation
- Induces a prominent interferon-driven transcriptional program
- Does not produce dominant terminal exhaustion signatures

These findings support a model of functional remodeling within pre-existing T cell populations during acute viral infection.

---

### Implementation

The project was independently implemented in:

**Python — Scanpy**
- Full scRNA-seq preprocessing and clustering
- T cell pseudobulk strategy
- DESeq2 integration

**R — Seurat**
- Reimplementation of the complete pipeline
- Harmony batch correction
- Refined CD4 and CD8 sub-analyses
- Cross-condition compositional comparisons

Both pipelines converge on consistent biological conclusions.

---

```
scRNA_Influenza/
│
├── Scanpy/
│   ├── scripts/
│   └── figures/
│
├── Seurat/
│   ├── scripts/
│   └── figures/
│
└── README.md

```
Large intermediate files (.h5ad, raw matrices) are excluded from version control.

---

### Key Take-Home Message

Acute influenza infection appears to reshape transcriptional programs within existing T cell populations, promoting interferon-driven activation and effector differentiation without evidence of dominant terminal exhaustion.