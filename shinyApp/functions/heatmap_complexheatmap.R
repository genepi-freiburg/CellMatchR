library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dplyr)

# function to flexibly adapt height of heatmap to number of genes
calculate_heatmap_height <- function(n_genes) {
  base_height <- 200
  height_per_gene <- 18
  min_height <- 400
  max_height <- 2000
  calculated <- base_height + (n_genes * height_per_gene)
  return(max(min_height, min(max_height, calculated)))
}

# heatmap function
create_complex_heatmap_shiny <- function(reference, sample, markers, cor_order, genesToPlot = NULL, 
                                         results_cor, selected_option, scale = FALSE,
                                         gene_celltype = NULL) {  # NEW parameter
  
  # Get filtered data
  if(selected_option == "counts"){
    data_filt <- reduce_data_counts(reference, sample)
  } else if(selected_option == "CPM") {
    data_filt <- reduce_data(reference, sample)
  }
  
  cellnames <- cor_order$celltypes
  
  # Filter by genes
  if(!is.null(genesToPlot) && length(genesToPlot) > 0) {
    print(paste("User selected", length(genesToPlot), "genes"))
    data_filt$ref_filt <- data_filt$ref_filt[rownames(data_filt$ref_filt) %in% genesToPlot,, drop = F]
    data_filt$test_filt <- data_filt$test_filt[rownames(data_filt$test_filt) %in% genesToPlot,, drop = F]
  } else {
    marker_gene_list <- markers$Gene.names
    # Limit to 100 for readability
    if(length(marker_gene_list) > 100) {
      marker_gene_list <- marker_gene_list[1:100]
      print("Limiting default display to 100 marker genes")
    }
    print(paste("Using default markers:", length(marker_gene_list), "genes"))
    data_filt$ref_filt <- data_filt$ref_filt[rownames(data_filt$ref_filt) %in% marker_gene_list,, drop = F]
    data_filt$test_filt <- data_filt$test_filt[rownames(data_filt$test_filt) %in% marker_gene_list,, drop = F]
  }
  
  print(paste("After filtering - ref rows:", nrow(data_filt$ref_filt), "test rows:", nrow(data_filt$test_filt)))
  
  # Rename columns
  cor_names <- as.data.frame(rownames(results_cor))
  ref_names <- cor_names[1:ncol(data_filt$ref_filt), ]
  test_names <- cor_names[ncol(data_filt$ref_filt) + 1 : ncol(data_filt$test_filt),]
  
  ref_filt <- data_filt$ref_filt
  colnames(ref_filt) <- paste0("ref_", ref_names)
  
  test_filt <- data_filt$test_filt
  colnames(test_filt) <- paste0("test_", test_names)
  
  # Identify which cellnames belong to test vs ref AFTER renaming
  test_col_order <- cellnames[cellnames %in% colnames(test_filt)]
  ref_col_order  <- cellnames[cellnames %in% colnames(ref_filt)]
  
  print(paste("test_col_order:", paste(test_col_order, collapse = ", ")))
  print(paste("ref_col_order:", paste(ref_col_order, collapse = ", ")))
  
  # Combine: sample (test) on LEFT, reference on RIGHT
  data <- cbind(
    test_filt[, test_col_order, drop = FALSE],
    ref_filt[, ref_col_order, drop = FALSE]
  )
  
  #============================================#
  # CELL TYPE GROUPING
  #============================================#
  
  if(!is.null(gene_celltype)) {
    present_genes <- rownames(data)
    
    # Filter cell type table to present genes
    gene_ct_filt <- gene_celltype %>%
      filter(Gene.names %in% present_genes)
    
    # Build expanded matrix: duplicate rows for genes in multiple cell types
    expanded_rows <- list()
    gene_labels <- c()
    cell_type_labels <- c()
    
    for(gene in present_genes) {
      cts <- gene_ct_filt$abbreviation[gene_ct_filt$Gene.names == gene]
      
      if(length(cts) == 0) {
        # Gene not in cell type table
        expanded_rows[[length(expanded_rows) + 1]] <- data[gene, , drop = FALSE]
        gene_labels <- c(gene_labels, gene)
        cell_type_labels <- c(cell_type_labels, "Other")
      } else if(length(cts) == 1) {
        expanded_rows[[length(expanded_rows) + 1]] <- data[gene, , drop = FALSE]
        gene_labels <- c(gene_labels, gene)
        cell_type_labels <- c(cell_type_labels, cts)
      } else {
        # Gene in multiple cell types - add .1, .2 etc suffix
        for(j in seq_along(cts)) {
          expanded_rows[[length(expanded_rows) + 1]] <- data[gene, , drop = FALSE]
          gene_labels <- c(gene_labels, paste0(gene, ".", j))
          cell_type_labels <- c(cell_type_labels, cts[j])
        }
      }
    }
    
    # Combine into matrix and sort by cell type
    data_expanded <- do.call(rbind, expanded_rows)
    rownames(data_expanded) <- gene_labels
    
    sort_order <- order(cell_type_labels)
    data_expanded <- data_expanded[sort_order, ]
    cell_type_labels <- cell_type_labels[sort_order]
    
    data <- data_expanded
    row_split_vector <- factor(cell_type_labels, levels = unique(cell_type_labels))
    
  } else {
    # No cell type info: cluster by ward.D2
    if(nrow(data) > 1) {
      ref_clust_genes <- hclust(dist(log2(ref_filt[,] + 1)), method = 'ward.D2')
      data <- data[ref_clust_genes$order,]
    }
    row_split_vector <- NULL
  }
  
  #============================================#
  # TRANSFORMATIONS
  #============================================#
  
  # Log2 transformation
  data_log2 <- log2(data + 1)
  
  # Scaling
  if(scale == TRUE) {
    normalize_expression_by_row <- function(x) {
      if (max(x) == 0) {
        return(rep(0, length(x)))
      } else {
        return((x - min(x)) / (max(x) - min(x)))
      }
    }
    scale_data <- as.data.frame(t(apply(data_log2, 1, normalize_expression_by_row)))
    colnames(scale_data) <- colnames(data_log2)
    data_log2 <- scale_data
  }
  
  # Convert to matrix
  expression_matrix <- as.matrix(data_log2)
  
  # Column groups: Sample FIRST (left), Reference SECOND (right)
  n_test <- length(test_col_order)
  n_ref  <- length(ref_col_order)
  column_groups <- factor(
    c(rep("Sample", n_test), rep("Reference", n_ref)),
    levels = c("Sample", "Reference")
  )
  
  # Color function
  col_fun <- colorRamp2(
    c(min(expression_matrix, na.rm = TRUE), 
      median(expression_matrix, na.rm = TRUE),
      max(expression_matrix, na.rm = TRUE)), 
    plasma(3)
  )
  
  # Adjust font size based on number of genes
  n_genes <- nrow(expression_matrix)
  row_fontsize <- if(n_genes > 100) 6 else if(n_genes > 50) 9 else 14
  print(paste("Creating heatmap with", n_genes, "genes, fontsize:", row_fontsize))
  
  #============================================#
  # CREATE HEATMAP
  #============================================#
  
  heatmap <- Heatmap(
    expression_matrix,
    name = ifelse(scale, "Scaled\nExpression", "log2(CPM + 1)"),
    col = col_fun,
    
    # Column settings
    column_split = column_groups,
    column_gap = unit(2, "mm"),
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    column_title = NULL,
    
    # Row settings - gene names on left
    row_names_side = "left",
    row_names_gp = gpar(fontsize = row_fontsize),
    show_row_names = TRUE,
    
    # Row grouping by cell type
    row_split = row_split_vector,
    row_gap = unit(1.5, "mm"),
    row_title = NULL,
    
    # Clustering - off when grouped by cell type
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    cluster_row_slices = FALSE,
    
    # Visual
    rect_gp = gpar(col = "white", lwd = 0.5),
    
    # Left annotation block
    left_annotation = if(!is.null(row_split_vector)) {
      rowAnnotation(
        "Cell type" = anno_block(
          gp = gpar(fill = "white", col = "black", lwd = 1),
          labels = levels(row_split_vector),
          labels_gp = gpar(col = "black", fontsize = 9, fontface = "bold"),
          labels_rot = 0
        )
      )
    } else NULL,
    
    # Top annotation
    top_annotation = HeatmapAnnotation(
      Group = column_groups,
      col = list(Group = c("Sample" = "#D55E00", "Reference" = "#0072B2")),
      annotation_name_gp = gpar(fontsize = 10),
      show_legend = TRUE
    ),
    
    # Legend
    heatmap_legend_param = list(
      title_position = "topleft",
      legend_height = unit(4, "cm"),
      legend_direction = "vertical"
    )
  )
  
  # Download matrix
  gene_names <- rownames(data_log2)
  download_matrix <- cbind(gene_names, data_log2)
  
  return(list(
    heatmap = heatmap, 
    download_matrix = download_matrix,
    n_genes = n_genes
  ))
}