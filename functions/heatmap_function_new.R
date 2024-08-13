library(plotly)


## data normalization between 0 and 1   
normalize_expression_by_row <- function(x) {
  if (max(x) == 0) {
    row <- rep(0, length(x)) #verbessern?
    return(row)
  } else {
    row <- (x-min(x))/(max(x)-min(x))
    return(row)
  }
}

heatmap_fct <- function(reference, sample, markers, cor_order, genesToPlot = NULL, results_cor, selected_option, scale = FALSE) {
  
  markers <- read.csv("data/Lake_2023_tubular_celltype_marker_genes.csv")
  
  if(selected_option == "counts"){
    data_filt <- reduce_data_counts(reference, sample)
  } else if(selection_option == "CPM") {
    data_filt <- reduce_data(reference, sample)
  }
  
  cellnames <- cor_order$celltypes               #vector containing celltypes in Spearman's cor order
  
  if(!is.null(genesToPlot)){ # use preset marker gene list in case no genes were selected 
    data_filt$ref_filt <-data_filt$ref_filt[rownames(data_filt$ref_filt) %in% genesToPlot,, drop = F]
    data_filt$test_filt <-data_filt$test_filt[rownames(data_filt$test_filt) %in% genesToPlot,, drop = F]
  } else { # if user has selected genes
    data_filt$ref_filt <-data_filt$ref_filt[rownames(data_filt$ref_filt) %in% markers$Gene.names,, drop = F]
    data_filt$test_filt <-data_filt$test_filt[rownames(data_filt$test_filt) %in% markers$Gene.names,, drop = F]
  }
  
  # problem that occured: when sample and reference have identical cellnames -> colnames of data_filt no longer identical to colnames of 
  # cellnames as cellnames derived from results_cor introduces "." to duplicate
  # solution: take rownames from results_cor which is already modified with "." if duplicates are present
  # and rename ref_filt and test_filt column according to rownames
  
  cor_names <- as.data.frame(rownames(results_cor))
  ref_names <- cor_names[1:ncol(data_filt$ref_filt), ]
  test_names <- cor_names[ncol(data_filt$ref_filt) + 1 : ncol(data_filt$test_filt),]
  
  ref_filt <- data_filt$ref_filt
  colnames(ref_filt) <- paste0("ref_", ref_names)  # add ref_ in front of ref_filt
  
  test_filt <- data_filt$test_filt
  colnames(test_filt) <- paste0("test_", test_names)

  data <- cbind(test_filt,ref_filt)
  
  
  # cluster and order genes (across reference cell types only)
  # Ward’s is said to be the most suitable method for quantitative variables.
  # if genesToPlot contains only 1 item, genes are NOT clustered
  # x-axis is ordered according to Spearman results order, y-axis according to cluster
  
  if(is.null(genesToPlot)) {
   ref_clust_genes <- hclust(dist(log2(ref_filt[,] + 1)), method = 'ward.D2')
   data <- data[, match(cellnames, colnames(data))]
   data <- data[ref_clust_genes$order,]
  } else if(length(genesToPlot) > 1) {
    ref_clust_genes <- hclust(dist(log2(ref_filt[,] + 1)), method = 'ward.D2')
    data <- data[, match(cellnames, colnames(data))]
    data <- data[ref_clust_genes$order,]
  } else {
    data <- data[, match(cellnames, colnames(data))]
   }

  # log2 transformation
  data_log2 <- log2(data + 1)
  
  if(scale == FALSE) {
    data_log2 <- data_log2
  } else if(scale == TRUE) {
    scale_data <- as.data.frame(t(apply(data_log2, 1, normalize_expression_by_row)))
    colnames(scale_data) <- colnames(data_log2)
    data_log2 <- scale_data
  }
  
  if(is.null(genesToPlot)){
    data_log2 <- data_log2
    
  } else {
    data_log2$Gene.names <- rownames(data_log2)
    data_select <- data_log2 %>%
      filter(Gene.names %in% genesToPlot)
    print(data_log2)
    data_select$Gene.names <- NULL
    data_log2 <- data_select
  }
  
  ## structure data_log2 for download
  gene_names <- rownames(data_log2)
  download_matrix <- cbind(gene_names, data_log2)
  
  map <- plotly::plot_ly(
    z = as.matrix(data_log2),
    x = colnames(data_log2),
    y = rownames(data_log2),
    type = "heatmap",
    zauto = TRUE
  )%>%
    plotly::layout(
      title = "",
      
      xaxis = list(
        title = "",
        categoryorder = "array", categoryarray = cor_order,
        mirror = TRUE,         
        autorange = TRUE,
        showline = FALSE,
        type = "category"
      ),
      yaxis = list(
        autorange = TRUE,
        hoverformat = ".2f",
        mirror = TRUE,
        showline = FALSE
      ),
      dragmode =  "select",
      hovermode = "compare"
    ) %>%
    plotly::config(
      toImageButtonOptions = list(format = "png",
                                  filename = "heatmap",
                                  height = 960,
                                  width = 1280,
                                  scale = 2)
    )
  
  map <- ggplotly(map)  %>% # here it seems layout can be formatted like in ggplot.
    layout(
      xaxis = list(tickfont = list(size = 10)), 
      yaxis = list(tickfont = list(size = 7)))
  
  return(list(map = map, download_matrix = download_matrix)) #in Shiny, only 1 object can be returned -> list
}


