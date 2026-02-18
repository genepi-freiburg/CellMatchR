# CellMatchRLive

**Cell type matching for kidney RNA-seq data - runs entirely in your browser**

**Live app:** https://genepi-freiburg.github.io/CellMatchR/

## About

CellMatchR compares gene expression profiles from bulk RNA-sequencing of kidney cell lines, primary cells, or whole tissue to single-cell RNA-seq references using Spearman correlation.

**Key features:**
- Runs entirely in browser using WebAssembly (no data uploaded to any server)
- Spearman correlation matching
- Heatmap with cell-type grouping of 40 curated tubular marker genes
- 4 reference datasets (mouse and human)
- Upload your own samples, references, or marker genes

## Technical Details

Built with [Shinylive](https://posit-dev.github.io/shinylive/) - R Shiny apps compiled to WebAssembly via [webR](https://docs.r-wasm.org/webr/latest/).

First load takes 30-60 seconds to download the R runtime (~30MB). Subsequent visits are faster due to browser caching.

The `docs/` folder contains the Shinylive export and is deployed to GitHub Pages via the workflow in `.github/workflows/deploy.yml`.
