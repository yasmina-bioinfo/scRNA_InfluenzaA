## 01_Setup_Seurat_Pipeline 
Packages required for Seurat scRNA-seq pipeline were installed manually prior to running the analysis scripts 

## 02. Create Seurat Objects
- Imported raw count matrices.
- Created individual Seurat objects per sample (IC excluded).

Output:
- seu_list_raw.rds

## 03. Quality Control Metrics
Computed standard QC metrics per sample, including total counts, detected features, and mitochondrial content.

Generated diagnostic plots to guide threshold selection:
- QC_scatter (nCount vs percent.mt)
- QC_violin per sample
- QC summaries were exported for transparency and reproducibility.

Thresholds applied (defined in Step 04):
- nFeature_RNA > 200
- nFeature_RNA < 4000
- percent.mt < 10

Outputs:
- Figures saved in Results/Figures/
- qc_summary_by_sample.csv saved in Results/qc/

## 04. Quality Filtering
Per-sample quality control was applied using the following thresholds:
- nFeature_RNA > 200 &
- nFeature_RNA < 4000 &
- percent.mt < 10
Cell loss remained below ~15% per sample, preserving dataset integrity while removing low-quality cells.

Outputs: 
- seu_list_filtered.rds and filtering_summary.csv
- Max 15% of cells lost 

## 05. Normalization 
Counts were log-normalized per sample to correct for differences in sequencing depth, enabling comparable gene expression analysis across cells.

Output
- seu_list_norm_hvg.rds

## 06. Merge & Harmony Batch Correction
Filtered samples were merged into a single Seurat object. After normalization and PCA, batch effects were corrected using Harmony, with sample identity specified as the grouping variable.

UMAP and clustering were subsequently performed on the Harmony-corrected embeddings to ensure that downstream structure reflects biological variation rather than technical effects.

## 07. Harmony Integration – PCA & UMAP Visualization

- Visualized cluster structure and sample distribution after Harmony integration.
- Assessed batch correction and biological signal preservation.

Outputs
- PCA_clusters.png
- PCA_samples.png
- UMAP_clusters.png
- UMAP_samples.png

## 08. Cluster Marker Identification 
Identify the main immune cell populations present in the dataset.
- Performed FindAllMarkers on seu_harmony_merged.rds (Harmony-integrated object).
- Identified differentially expressed genes for each cluster.
- Selected the top 20 markers per cluster based on adjusted p-value and logFC.
- Used canonical immune signatures to support downstream cell-type annotation.

Output
- markers_all_clusters.tsv

## 09. Broad CellType Annotation + Global UMAP

- Assigned family-level Cell_Annotation labels to all clusters (T/NK, B/Plasma, Myeloid, Neutrophils, Others).
- Validated full cluster coverage and exported a summary table.
- Generated a global UMAP colored by lineage using a fixed palette for consistency across analyses.

Outputs
- cell_annotation_summary.tsv
- UMAP_global_by_Cell_Annotation.png
- Updated Seurat object with Cell_Annotation

## 10 & 11. T cell Subsetting and Annotation (from seu_tcells_annotated.rds)
- Sub-clustering of the global T cell compartment
- Assignment of lineage and functional identities based on canonical markers

Outputs
- DotPlot_Tcells_Global.png
- UMAP_tcells_by_cluster.png

Cluster annotation logic
- CD4_Naive_TCM: high IL7R, CCR7, LEF1 → naïve / central memory phenotype
- CD8_Effector / CD8_Terminal_Cytotoxic: NKG7, GNLY, GZMB, PRF1, - - FGFBP2 = cytotoxic differentiation
- ISG_high_T: ISG15, IFIT2, IFI27, OAS1, OASL → interferon-responsive subset
- GammaDelta_T: TRDC, TRGC1, TRGC2 → γδ T cells
MAIT_like: KLRB1 with TRAV1-2 / SLC4A10 → MAIT-like population

## 12 Focus on CD8 T cells
To characterize functional CD8 states following Influenza A infection, we analyzed the 12 clusters obtained from the T cell subset.
CD8 T cells display a functional continuum ranging from memory-like states (CCR7, IL7R) to highly cytotoxic effectors (GZMB, PRF1, NKG7, GNLY). An interferon-responsive subset (ISG15, STAT1) is also observed, reflecting antiviral activation.
Classical exhaustion markers (PDCD1, LAG3, TOX) are detectable but remain low to moderate and do not define a dominant cluster, consistent with acute viral infection rather than chronic T cell dysfunction.
These results indicate a robust and diversified antiviral CD8 response, characterized by active effector differentiation without evidence of terminal exhaustion.

Saved as:
- DotPlot_CD8_by_cluster.png
- UMAP_CD8_annotated.png

## 13 Focus on CD4 T cells 
To characterize functional CD4 states following Influenza A infection, we analyzed the clusters obtained from the CD4 T cell subset.
CD4 T cells display a structured continuum ranging from naive/memory-like states (CCR7, IL7R, TCF7) to Th1-like populations (CXCR3, IFNG), as well as a smaller cytotoxic-like subset expressing GZMB and PRF1. An interferon-high population (ISG15, STAT1) is clearly segregated, reflecting strong antiviral activation.

FOXP3 and CTLA4 expression identify a discrete regulatory T cell (Treg) cluster. No dominant exhaustion program is observed, consistent with an acute antiviral immune response rather than chronic dysfunction.

These results indicate a diversified and coordinated CD4 response during Influenza A infection, combining helper polarization, cytotoxic potential, regulatory balance, and interferon-driven activation.

Saved as:
- DotPlot_CD4_by_cluster.png
- UMAP_CD4_annotated.png

**Integrated CD4–CD8 perspective**

Together, CD4 and CD8 T cell analyses reveal a coordinated antiviral immune response following Influenza A infection.

CD8 T cells predominantly adopt cytotoxic effector programs, reflecting direct antiviral activity. In contrast, CD4 T cells display a broader functional diversification, including Th1 polarization, regulatory T cells, cytotoxic-like subsets, and a prominent interferon-responsive population.

Importantly, neither compartment shows evidence of dominant terminal exhaustion, supporting the interpretation of an acute, active immune response rather than chronic dysfunction.

These findings highlight complementary roles of CD4 and CD8 T cells in orchestrating effective antiviral immunity.

## 14 CD4_Comparison by Condition (HC,IA, PHC, PI)
To assess how CD4 T cell states vary across clinical groups (HC, IA, PHC, PI), we examined compositional shifts and transcriptional signatures.

Healthy controls (HC) are dominated by naïve/memory-like CD4 cells (CCR7, IL7R, TCF7), reflecting a basal immune state. In contrast, Influenza A infection (IA) induces a marked expansion of IFN-high CD4 cells, characterized by strong ISG15 and STAT1 expression. A similar interferon-driven signature is observed in pregnant infected individuals (PI). Pregnant healthy controls (PHC) display an intermediate profile, retaining a substantial naïve/memory compartment with moderate Th1-like representation.

Barplot analysis confirms the dominance of IFN-high CD4 states in infected groups (IA and PI), whereas HC and PHC maintain greater subset diversity. DotPlot and FeaturePlot analyses validate these compositional changes at the gene level.

Overall, Influenza A infection reshapes the CD4 compartment toward a pronounced interferon-associated activation state, contrasting with the naïve-dominant baseline observed in healthy individuals.

Saved as:
- UMAP_CD4_HC_IA_PHC_PI.png
- Barplot_CD4_states_HC_IA_PHC_PI.png
- DotPlot_CD4_HC_IA_PHC_PI.png
- FeaturePlot_CD4_global.png
- FeaturePlot_CD4_split_by_group.png

## 15 CD8_Comparison by Condition (HC, IA, PHC, PI)
To assess how CD8 functional states vary across conditions, we quantified their distribution and gene expression profiles in HC, IA, PHC and PI samples using UMAP visualization, proportional barplots, DotPlots, and FeaturePlots.

Across infected conditions (IA and PI), CD8 T cells show a strong enrichment of cytotoxic programs, characterized by increased expression of GZMB, PRF1, NKG7, and GNLY, consistent with active antiviral effector differentiation. Interferon-responsive genes (ISG15, STAT1) are also elevated, reflecting acute immune activation.

In contrast, healthy conditions (HC and PHC) display a higher proportion of early memory/transitional CD8 states, with stronger expression of IL7R, CCR7, and TCF7, indicative of a basal or pre-effector phenotype.

Exhaustion-associated markers (PDCD1, LAG3, TOX) remain detectable but do not define a dominant dysfunctional population, supporting the interpretation of an acute, functional antiviral response rather than chronic exhaustion.

Overall, these results confirm a condition-dependent shift from memory-like CD8 states in healthy individuals toward cytotoxic effector dominance in infected groups.
Saved as:
- UMAP_CD8_HC_IA_PHC_PI.png
- Barplot_CD8_states_HC_IA_PHC_PI.png
- DotPlot_CD8_HC_IA_PHC_PI.png
- FeaturePlot_CD8_global.png
- FeaturePlot_CD8_split_by_group.png

## Fuctional Remodeling of T cCells Across Conditions 

| Condition | CD4 State Profile                             | CD8 State Profile                      | Putative Functional Signal                    |
| --------- | --------------------------------------------- | -------------------------------------- | --------------------------------------------- |
| **HC**    | Predominantly Naive / Memory                  | Early Memory / Transitional            | Baseline immune surveillance                  |
| **PHC**   | Naive / Memory with mild Th1 tendency         | Memory / Transitional                  | Immune preparedness without strong activation |
| **IA**    | Increased Th1-like and IFN-responsive subsets | Expansion of Cytotoxic Effector state  | Active antiviral response                     |
| **PI**    | Sustained Th1-like / IFN-high signature       | Dominant Cytotoxic Effector population | Ongoing effector differentiation              |

## Cross-Lineage Interpretation

| Axis                     | CD4                                 | CD8                                 |
| ------------------------ | ----------------------------------- | ----------------------------------- |
| Effector differentiation | Appears induced during infection    | Strongly expanded during infection  |
| IFN response             | Distinct IFN-high subset detectable | Integrated within cytotoxic program |
| Exhaustion markers       | Detectable but not dominant         | Detectable but not dominant         |
| Functional plasticity    | Likely preserved                    | Likely preserved                    |

## Take-home interpretation

**These observations suggest that acute Influenza infection may promote coordinated CD4 Th1 polarization and CD8 cytotoxic differentiation, without clear evidence of terminal exhaustion.**