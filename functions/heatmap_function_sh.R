  #    reference <- read.csv("../data/Ransick_CPM_Gene.names.csv")
  # markers <- read.csv("../data/Lake_2023_tubular_celltype_marker_genes.csv")
  #  rownames(reference) <- reference$Gene.names
  #  chen <- read.csv("../../../../05_Output/01_Counts_CPM_TPM/01_Mona_processed_tables/Cpm_Chen_median.csv")
  #  
  #  sample <- chen[, c(1, 15)]
  #  sample <- read.csv(".../../../testing/3D_WT.csv")
  # source("reduce_data.R")
  #  source("cor_data.R")

#library(heatmaply)
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
   
  
   
heatmap_fct <- function(reference, sample, markers) {
  
  data_filt <- reduce_data(reference, sample, markers)
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
  plot(ref_clust_genes)

  #ref_clust_cells <- hclust(dist(log2(t(ref_filt[ref_clust_genes$order,] + 1))), method = 'average')
  
  # order genes in data according to clusters
  data <- data[ref_clust_genes$order,]

  # log2 transformation
  data_log2 <- log2(data + 1)
  
  #scale by row
  scale_log2 <- t(apply(data_log2, 1, normalize_expression_by_row))

  map <- plotly::plot_ly(
     z = as.matrix(scale_log2),
     x = colnames(data_log2),
     y = rownames(data_log2),
    
    #height = expression_heatmap_height(),
    # xgap = 4,
    # ygap = 4,
    type = "heatmap",
    #colorscale= color_scale,
    zauto = TRUE,
    #zmin = colorscale_min,
    #zmax = colorscale_max,
    #reversescale = FALSE,
    # colorbar = list(
    #   title = "Expression",
    #   lenmode = "pixels",
    #   len = 120
    # )
  )%>%
    plotly::layout(
      title = "",
      
      xaxis = list(
        title = "",
        mirror = TRUE,         # mirrors axis lines to opposite plto site
        autorange = TRUE,
        #dtick = 1,
        showline = FALSE,
        type = "category"
        # categoryorder = "array",    # "category ascending" does not work for heatmaps
        # categoryarray = sort(unique(expression_heatmap_data()$clusters))
      ),
      yaxis = list(
        #title = "Gene(s)",
        autorange = TRUE,
        hoverformat = ".2f",
        mirror = TRUE,
        #dtick = 1,
        showline = FALSE
      ),
      dragmode =  "select",
      hovermode = "compare"
    )
  
  map
  map <- ggplotly(map)  %>% # here it seems layout can be formatted like in ggplot.
    layout(
      xaxis = list(tickfont = list(size = 10)), 
      yaxis = list(tickfont = list(size = 7)))
  map

  
  # heatmap
   
    # map <- plot_ly(
    #    z = as.matrix(scale_log2),
    #    x = colnames(scale_log2),
    #    y = rownames(scale_log2),
    #   # height = 600,
    #   # width = 600,
    #    type = "heatmap",
    #   yaxis_nticks=dim(data_filt)[1]
    #  )
    # 
    # map <- ggplotly(map)
  
  return(map)
}



##----------------------------------------------------------------------------##
## DATA: EXPRESSION HEATMAP
##----------------------------------------------------------------------------##




##----------------------------------------------------------------------------##
## PLOT: EXPRESSION HEATMAP
##----------------------------------------------------------------------------##
# 
# 
# # function for normalizing gene expression averages/per cluster for each gene (considering that some genes might have zero expression for all clusters)
# normalize_expression_by_row <- function(x) {
#   if (max(x) == 0) {
#     row <- rep(0, length(x))
#     return(row)
#   } else {
#     row <- (x-min(x))/(max(x)-min(x))
#     return(row)
#   }
# }
# 
# expression_heatmap_data <- reactive({
#   req(
#     #input[["global_cluster_select"]],
#     #input[["expression_heatmap_rescaling_select"]]
#   )
# 
#   clusters_to_display <- input[["global_cluster_select"]]
#   rescale <- input[["expression_heatmap_rescaling_select"]]
# 
#   # find which genes to display
#   reduce_data(reference(), sample(), geneset())
# 
#   if (type_exists() & sample_exists()) {
#     cells_to_display <-
#       which(sample_data()$cells$cluster %in% clusters_to_display() &
#               sample_data()$cells$type %in% types_to_display() &
#               sample_data()$cells$sample %in% samples_to_display())
#   } else if (type_exists()) {
#     cells_to_display <-
#       which(sample_data()$cells$cluster %in% clusters_to_display() &
#               sample_data()$cells$type %in% types_to_display())
#   } else if (sample_exists()) {
#     cells_to_display <-
#       which(sample_data()$cells$cluster %in% clusters_to_display() &
#               sample_data()$cells$sample %in% samples_to_display())
#   } else {
#     cells_to_display <-
#       which(sample_data()$cells$cluster %in% clusters_to_display())
#   }
# 
#   # in case wrong cluster names are loaded in control:
#   if (is.null(cells_to_display)) {
#     clusters_to_display <- sample_data()$cluster_names
#     cells_to_display <- which(
#       sample_data()$cells$cluster %in% clusters_to_display)
#   }
# 
#   cluster_matrix <- sample_data()$cells[cells_to_display, c("cluster", "cell_barcode")]
#   genes <- genesToPlot()$genes_to_display_present
#   expression_matrix <- t(sample_data()$expression[genes, cells_to_display, drop=FALSE]) # drop false ==> keep dcgMatrix format also for one col
#   # convert from sparse format
#   expression_matrix <- as.matrix(expression_matrix)
#   hmap_dataframe <- cbind(cluster_matrix, expression_matrix)
#   hmap_dataframe[2] <- NULL  # remove cell barcodes
# 
#   if( ncol(hmap_dataframe) > 1) {
# 
#     # sort by cluster (code for cell informaton currently commented out)
#     hmap_dataframe <- hmap_dataframe[order(hmap_dataframe$cluster),]
# 
#     # dataframe for expression per cell
#     hmap_dataframe_cells <- t(hmap_dataframe)
#     colnames(hmap_dataframe_cells) <- hmap_dataframe_cells[1,]
#     hmap_dataframe_cells <- hmap_dataframe_cells[-1,,drop=F]
#     hmap_matrix_cells <- data.matrix(hmap_dataframe_cells)
#     class(hmap_matrix_cells) <- "numeric"
# 
#     # calculate average expression per cluster
#     hmap_dataframe <- hmap_dataframe %>% group_by(cluster) %>% summarise_all(funs(mean))
# 
#     # remove cluster column
#     hmap_matrix <- t(as.matrix(hmap_dataframe[-1]))
# 
#     #normalize if selected and more than one cluster selected
#     if(rescale == "Rescale per gene" && length(clusters_to_display) > 1) {
#       hmap_matrix <- t(apply(hmap_matrix, 1, normalize_expression_by_row))
#       #  hmap_matrix <- t(apply(hmap_matrix, 1, function(x){(x-min(x))/(max(x)-min(x))}))
#     }
# 
#   } else {
#     hmap_matrix <- NULL
#     hmap_matrix_cells <- NULL
#   }
# 
#   # generate ouput data
#   hmap_data <- list()
#   hmap_data[["clusters"]] <- as.vector(hmap_dataframe$cluster)
#   hmap_data[["genes"]] <- as.vector(colnames(hmap_dataframe))
#   hmap_data[["genes"]] <- hmap_data[["genes"]][-1]    # remove "cluster"
#   hmap_data[["av_expression"]] <- hmap_matrix
#   hmap_data[["expression_cells"]] <- hmap_matrix_cells
# 
#   return(hmap_data)
# })
