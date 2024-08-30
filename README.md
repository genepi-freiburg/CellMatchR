# CellMatchR
# CellMatchR is an efficient tool to statistically compare gene expression profiles derived from  # RNA sequencing of whole tissues, primary cells and cell lines of the kidney.
# For references, we are using published kidney single cell RNA sequencing datasets. For 
# statistical comparisons,
# various approaches were tested for best performance and we settled for simple Spearman’s 
# correlation coefficients and Euclidean distance.
# Matching can be run across all genes in the dataset or across a chosen set of genes, such as 
# tubular kidney marker genes.
# Reference datasets were derived from the published manuscript. 
# Pseudobulk counts of each cell type and gene were created by adding up the counts for each gene 
# across all single cells annotated to the same cell. 
# Counts data were then normalized by calculating the counts per million (CPM) of each pseudobulk
# cell type.
