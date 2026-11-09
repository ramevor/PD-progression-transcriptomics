# A Transcriptomic Signature of Parkinson's Disease Progression

R code for the analyses in *"A Transcriptomic Signature of Parkinson's Disease Progression: Evidence from the PPMI Study"* (Amevor et al.).

This repository contains the full analysis pipeline: cohort construction, univariate Cox gene screening, penalized Cox regression (LASSO / MCP / SCAD), gene set enrichment analysis (GSEA), model evaluation (cross-validated C-index and Integrated Brier Score), and figure generation.

## Data availability

The Parkinson's Progression Markers Initiative (PPMI) data are **access-controlled and are not included in this repository**. To reproduce the analysis, register and download the data from https://www.ppmi-info.org.

You will need:
- `counts_matrix_TPM.csv` — whole blood RNA-seq TPM matrix
- `PPMI_Curated.csv` — curated clinical/demographic data
- MSigDB gene set files (`.gmt`) for KEGG, Reactome, and GO Biological Process, available from https://www.gsea-msigdb.org/gsea/msigdb

## Setup

Place the data files as follows (or edit the paths in the **USER CONFIG** block at the top of the script):

```
data/
  counts_matrix_TPM.csv
  PPMI_Curated.csv
  gmt/
    c2.cp.kegg.v2022.1.Hs.symbols.gmt
    c2.cp.reactome.v2026.1.Hs.symbols.gmt.txt
    c5.go.bp.v2026.1.Hs.symbols.gmt.txt
output/          # created automatically for figures and tables
```

## Required R packages

```r
install.packages(c("readr","dplyr","stringr","lubridate","survival","survminer",
                   "purrr","broom","tibble","glmnet","ncvreg","caret","ggplot2",
                   "tableone","gridExtra","openxlsx","pec","ggVennDiagram"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("org.Hs.eg.db","AnnotationDbi","fgsea"))
```

## Running

Open `PD_Analysis_GitHub.R`, set the paths in the **USER CONFIG** block, and run the script from top to bottom. Outputs (figures, tables, GSEA results) are written to the `output/` folder.

## Pipeline overview

| Section | Description |
|---------|-------------|
| 0 | Packages, theme, output setup |
| 1 | Data loading, TPM filtering, log2 transform |
| 2 | Cohort definition (HY 1–2 at baseline → HY ≥ 3 transition) |
| 3 | Descriptive statistics (Table 1), Kaplan–Meier curves |
| 4 | Univariate and age/sex/race-adjusted Cox gene screening |
| 5 | Penalized selection (LASSO / MCP / SCAD), Venn diagram |
| 6 | Cross-validated C-index for nested models |
| 7 | Forest plot of MCP gene effects |
| 8 | Risk-group stratification and Kaplan–Meier plot |
| 9 | Integrated Brier Score (calibration) |
| GSEA | Gene set enrichment across KEGG, Reactome, GO BP |

## Citation

If you use this code, please cite the associated paper (citation to be added upon publication).
