# CellMatchR

CellMatchR is an efficient tool to compare gene expression profiles derived from RNA sequencing of whole tissues, primary cells, and cell lines of the kidney using statistical and machine learning methods.

## References

We are using published kidney single cell RNA sequencing datasets as references. For comparisons, various approaches were tested for optimal performance, we settled on simple Spearman’s correlation coefficients and Euclidean distance for statisical comparisons and TabPFN as machine learning model.

## Matching

Matching can be performed across all genes in the dataset or across a chosen set of genes, such as tubular kidney marker genes. Reference datasets were derived from a published manuscript.

## Data Processing

Pseudobulk counts of each cell type and gene were created by summing the counts for each gene across all single cells annotated to the same cell type. The counts data were then normalized by calculating the counts per million (CPM) for each pseudobulk cell type.


