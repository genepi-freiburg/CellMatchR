#-----------------------------------------------#
#  Heatmap Code
#----------------------------------------------#

library(dplyr)


# create reactive list containing gene names of reduced data objects

red_data <- reactive({
  if(is.null(rv_sample$data) & input$demo == "--") {
    return(NULL)
  } else if(input$button > 0) {
    red_data <- reduce_data(reference(), sample(), geneset())
    ref_filt <- red_data$ref_filt
    gene_name_list <<- as.list(rownames(ref_filt))
    return(gene_name_list)
  }
  else{
    return(NULL)
  }
})

# updates selectizeInput to contain gene names

observe({
  updateSelectizeInput(session, "expression_genes_input",
                       choices = red_data(),
                       server = TRUE)
})

# when clicking reset button, selected genes are removed

observeEvent(input$reset_genes, {
  updateSelectizeInput(session, "expression_genes_input",
                       choices = red_data(),
                       selected = NULL,
                       server = TRUE)
  
})

# create reactive value list for genes chosen by user

selected_genes <- reactiveValues(data = NULL, gene_list = list())

observe({
  selected_genes$gene_list <- input$expression_genes_input
  })

genesToPlot <- reactive({
  selected_genes$gene_list 
})

# create reactive object containing Spearman results table ordered by descending rho, add ref_ and test_

cor_order <- reactive({
  ordered_cor <- results_cor() %>%
    mutate(celltypes = ifelse(datatype == "Median of sample(s) vs. reference", paste0("ref_", celltypes),
                              ifelse(datatype == "Median of sample(s) vs. sample(s)", paste0("test_", celltypes), celltypes))) %>%
    arrange(desc(rho))
  return(ordered_cor)
  })

# create heatmap output

results_map <- reactive({
  if(input$scale == "no scaling") {
    
    results_HM <- ({
      heatmap_fct(reference()[, c("Gene.names", input[["reference_celltype_select"]])],
                  sample()[, c("Gene.names", input[["test_celltype_select"]])],
                  markers,
                  cor_order(),
                  genesToPlot(),
                  results_cor(),      # results_cor provides colnames for heatmap dataframe 
                  selected_option(), #counts or cpm
                  scale = FALSE)
    })
   
  } else if(input$scale == "scale by row") {
    results_HM <- ({
      heatmap_fct(reference()[, c("Gene.names", input[["reference_celltype_select"]])],
                  sample()[, c("Gene.names", input[["test_celltype_select"]])],
                  markers,
                  cor_order(),
                  genesToPlot(),
                  results_cor(),
                  selected_option(),
                  scale = TRUE)
    })
  }
  return(results_HM)
})

output$heatmap <- renderPlotly({
  if(input$button > 0) {
    results_map()$map
  }
})
