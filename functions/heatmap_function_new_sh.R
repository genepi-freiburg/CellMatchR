library(plotly)


heatmap_fct <- function(reference, sample, markers, red_data = NULL) {
  
  markers <- read.csv("data/Lake_2023_tubular_celltype_marker_genes.csv")
  data_filt <- reduce_data(reference, sample)

  if(!is.null(red_data)){ # use preset marker gene list in case no genes were selected
    data_filt$ref_filt <-data_filt$ref_filt[rownames(data_filt$ref_filt) %in% red_data,]
    data_filt$test_filt <-data_filt$test_filt[rownames(data_filt$test_filt) %in% red_data,]
  } else { # if user has selected genes
    data_filt$ref_filt <-data_filt$ref_filt[rownames(data_filt$ref_filt) %in% markers$Gene.names,]
    data_filt$test_filt <-data_filt$test_filt[rownames(data_filt$test_filt) %in% markers$Gene.names,]
  }
    
  ref_filt <- data_filt$ref_filt
  
  # order cell types according to spearman  == x axis will always be ordered as matching results
  cor <- cor_data(reduce_data(reference, sample, markers))
  cor <- cor[1:dim(ref_filt)[2],]
  ref_filt <- ref_filt[,order(cor$rho, decreasing = T)]
  
  test_filt <- data_filt$test_filt
  
  data <- cbind(test_filt,ref_filt)
  
  # cluster and order genes (across reference cell types only)
  #  Ward’s is said to be the most suitable method for quantitative variables.
  ref_clust_genes <- hclust(dist(log2(ref_filt[,] + 1)), method = 'ward.D2')
  
  #ref_clust_cells <- hclust(dist(log2(t(ref_filt[ref_clust_genes$order,] + 1))), method = 'average')
  
  # order genes in data according to clusters
  data <- data[ref_clust_genes$order,]
  
  # log2 transformation
  data_log2 <- log2(data + 1)
  
  if(is.null(red_data)){
    data_log2 <- data_log2
    
  } else {
    data_log2$Gene.names <- rownames(data_log2)
    data_select <- data_log2 %>%
      filter(Gene.names %in% red_data)
    data_select$Gene.names <- NULL
    data_log2 <- data_select
  }
  
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
        mirror = TRUE,         # mirrors axis lines to opposite plto site
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
    )
  
  map <- ggplotly(map)  %>% # here it seems layout can be formatted like in ggplot.
    layout(
      xaxis = list(tickfont = list(size = 10)), 
      yaxis = list(tickfont = list(size = 7)))
  
  return(map)
}
