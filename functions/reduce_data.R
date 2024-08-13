# reduce_data, Shiny


library(openxlsx)
library(dplyr)
library("ggfortify")
library("ggplot2")
library("ggrepel")
library("ggpubr")
library("ggdendro")
library("gridExtra")


reduce_data <- function(ref, test, marker_genes = NULL){
 
  test$Gene.names <- toupper(test$Gene.names)
  ref$Gene.names <- toupper(ref$Gene.names)
  
  ref <- ref[!(duplicated(ref$Gene.names)),]
  test <- test[!(duplicated(test$Gene.names)),]
  
  test_filt <- na.omit(test)
  ref_filt <- na.omit(ref)
  
  if(is.null(marker_genes)) {
    gene_overlap <- data.frame(Gene.names = intersect(test_filt$Gene.names, ref_filt$Gene.names))
  } else {
    gene_overlap <- intersect(test_filt$Gene.names, ref_filt$Gene.names)
    marker_genes <- toupper(marker_genes$Gene.names)
    gene_overlap <- data.frame(Gene.names = intersect(gene_overlap, marker_genes))
  }
  
  test_filt <- left_join(gene_overlap, test, by = c("Gene.names"))
  ref_filt <- left_join(gene_overlap, ref,  by = c("Gene.names"))

  rownames(test_filt) <- test_filt$Gene.names
  test_filt$Gene.names <- NULL
  
  rownames(ref_filt) <- ref_filt$Gene.names
  ref_filt$Gene.names <- NULL
  
  test_median <- apply(X = test_filt, MARGIN = 1, FUN = median)
  
  QC = round((nrow(test_filt) / nrow(test)) * 100, 2)
  QC = paste0("(", QC, "% of ",  nrow(test), " genes used)")
  
  no_NA <- nrow(test) - nrow(test_filt)         #how many NAs were in the dataset
  
  output <- list(ref_filt = ref_filt,
                 test_filt = test_filt,
                 test_median = as.numeric(test_median),
                 QC = QC,
                 no_NA = no_NA)
 
  return(output)
}

reduce_data_counts <- function(ref, test, marker_genes = NULL){
  
  cpm_test <- as.data.frame(apply(test[, -grep("Gene.names", colnames(test))], MARGIN = 2, FUN = function(x) {x / sum(x) * 1000000}))
  cpm_test$Gene.names <- test$Gene.names
  test <- cpm_test
  
  test$Gene.names <- toupper(test$Gene.names)
  ref$Gene.names <- toupper(ref$Gene.names)
  
  ref <- ref[!(duplicated(ref$Gene.names)),]
  test <- test[!(duplicated(test$Gene.names)),]
  
  test_filt <- na.omit(test)
  ref_filt <- na.omit(ref)
  
  if(is.null(marker_genes)) {
    gene_overlap <- data.frame(Gene.names = intersect(test_filt$Gene.names, ref_filt$Gene.names))
  } else {
    gene_overlap <- intersect(test_filt$Gene.names, ref_filt$Gene.names)
    marker_genes <- toupper(marker_genes$Gene.names)
    gene_overlap <- data.frame(Gene.names = intersect(gene_overlap, marker_genes))
  }
  
  test_filt <- left_join(gene_overlap, test, by = c("Gene.names"))
  ref_filt <- left_join(gene_overlap, ref,  by = c("Gene.names"))
  
  rownames(test_filt) <- test_filt$Gene.names
  test_filt$Gene.names <- NULL
  
  rownames(ref_filt) <- ref_filt$Gene.names
  ref_filt$Gene.names <- NULL
  
  test_median <- apply(X = test_filt, MARGIN = 1, FUN = median)
  
  QC = round((nrow(test_filt) / nrow(test)) * 100, 2)
  QC = paste0("(", QC, "% of ",  nrow(test), " genes used)")
  
  no_NA <- nrow(test) - nrow(test_filt)         #how many NAs were in the dataset
  
  output <- list(ref_filt = ref_filt,
                 test_filt = test_filt,
                 test_median = as.numeric(test_median),
                 QC = QC,
                 no_NA = no_NA)
  
  return(output)
}