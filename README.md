# CellMatchR

## Overview

CellMatchR is a computational framework for matching bulk RNA-seq data - primarily from cell lines - to single-cell RNA-seq (scRNA-seq) reference datasets from the kidney. This repository accompanies the publication **Transcriptome-based cell type assignment for kidney cell culture models** from Schoberth, Boehm *et al.*. *(under review)*.

CellMatchR was generated to help researchers classify and annotate bulk RNA-seq samples from kidney cell lines by leveraging curated kidney scRNA-seq reference atlases, using both correlation-based approaches and machine learning (TabPFN).

---

## Contents

| Folder | Description |
|--------|-------------|
| `shinyApp/` | Interactive Shiny application for bulk-to-single-cell matching based on gene expression ranks |
| `TabPFN_scripts/` | User manual and scripts for TabPFN-based cell type classification |

---

## Shiny Application

The CellMatchR Shiny app allows users to upload bulk RNA-seq data and match it against curated kidney scRNA-seq references interactively — no coding required.

🔗 **[Launch CellMatchR App](https://epi.uniklinik-freiburg.de/cellmatchr)** — https://epi.uniklinik-freiburg.de/cellmatchr

### Features
- Upload your own bulk RNA-seq expression data
- Match against a curated kidney scRNA-seq reference atlas
- Method: Spearman's correlation
- Heatmap to visualize marker gene expression of samples and references
- Downloadable output tables

### Input Format
- comma-separated (.csv) or excel (.xlsx) expression matrix
- Rows: genes (HGNC symbols) in column "Gene.names", Columns: samples

---

## TabPFN User Manual

The `TabPFN/` folder contains a step-by-step guide to implement TabPFN for cell type matching of bulk RNA-seq data, including:

- Installation and environment setup
- Input data formatting requirements
- Running the TabPFN classifier
- Interpreting and exporting results

TabPFN is a transformer-based model for tabular data classification and provides a complementary machine learning approach to correlation-based matching.

---

## Requirements

### Shiny App (if installed locally)
- R (>= 4.0)
- See `shinyApp/` for required R packages

### TabPFN
- Python (>= 3.8)
- See `TabPFN/` user manual for full setup instructions

---

## Citation

If you use CellMatchR in your research, please cite:

> **Schoberth, Boehm et al.**. *Transcriptome-based cell type assignment for kidney cell culture models*. [Journal], [Year]. [DOI]

---

## Contact

For questions, bugs, or feature requests please open a [GitHub Issue](https://github.com/yourusername/CellMatchR/issues) or contact **[your email]**.


