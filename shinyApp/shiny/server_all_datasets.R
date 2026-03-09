#------------------------------------------------------------------------------#
# load all datasets, read updated datasets, update dataset selection list, results
#------------------------------------------------------------------------------#

reference <- reactive({
  switch(input$ref,
         "Ransick et al. (mouse, recommended)" = {
           result <- read.csv("data/Ransick_CPM_2025-08-14.csv")
         },
         "Park et al. (mouse)" = {
           result <- read.csv("data/Park_CPM_2025-08-14.csv")
         },
         "Lake et al. (human)" = {
           result <- read.csv("data/KPMP_CPM_2025-08-14.csv")
         },
         "Zhang et al. (human)" = {
          result <- read.csv("data/Zhang_CPM_2025-08-14.csv")
         }
  )
  if(is.null(result)) {
    message("Error reading the selected reference dataset.")
  }
  result
})

rv_sample <- reactiveValues(data = NULL)

observe({
  req(input$test)
  file_extension <- tools::file_ext(input$test$name)
  if (tolower(file_extension) == "csv") {
    rv_sample$data <- as.data.frame(read.csv(input$test$datapath))
  } else if (tolower(file_extension) == "xlsx") {
    rv_sample$data <- as.data.frame(readxl::read_excel(input$test$datapath))
  } else {
    stop("Unsupported file type. Please upload a CSV or Excel file.")
  }
})

observeEvent(input$test,{
  updateSelectInput(session, "demo", selected = "--")
})

observeEvent(input$reset_test, {
  rv_sample$data <- NULL
  reset(id = "test", asis = FALSE)
})

sample <- reactive({
  req(input$demo)

  if (!is.null(rv_sample$data) & input$demo == "--") {  
    sample <- rv_sample$data
  #   file_extension <- tools::file_ext(input$test$name)
  #   #updateSelectInput(session, "demo", selected = "--")
  #   
  #   if (tolower(file_extension) == "csv") {
  #     return(as.data.frame(read.csv(input$test$datapath)))
  #   } else if (tolower(file_extension) == "xlsx") {
  #     return(as.data.frame(readxl::read_excel(input$test$datapath)))
  #   } else {
  #     stop("Unsupported file type. Please upload a CSV or Excel file.")
  #   }
  #}
  # } else if(is.null(rv_sample$data) & input$demo == "--") {
  #   stop("Please select or upload test samples.")
  } else {
    selected_demo <- input$demo
    
    if (selected_demo == "kidney primary cells - select cell types below") {
      return(read.csv("data/CPM_Chen_primary_cells.csv"))
    } else if (selected_demo == "HK-2 proximal tubule cell line") {
      return(read.csv("data/TPM_HK2_khundmiri.csv"))
    } else if (selected_demo == "mIMCD-3 cell line") {
      return(read.csv("data/mIMCD_WT_10d_Westermann.csv"))
    }
  }
})


geneset <- reactive({
  selected_genes <- input$genes
  if(selected_genes == "all genes") {
    return(NULL)
  } else if(selected_genes == "kidney marker genes") {
    return(read.csv("data/marker_genes_all_2025-07-16.csv"))
  } else if(selected_genes == "tubular marker genes") {
    return(read.csv("data/marker_genes_tubular_2025-07-16.csv"))
  } else if (selected_genes %in% names(rv$gene_lists)) {
    return(rv$gene_lists[[selected_genes]])
  }
  })

# function to read uploaded genelist
read_gene_list <- function(genelist_input) {
  ext <- tools::file_ext(genelist_input$name)
  if(ext == "csv") {
    return(read.csv(genelist_input$datapath))
  } else if (ext == "xlsx") {
    return(readxl::read_excel(genelist_input$datapath))
  } else {
    return(NULL)
  }
} 
# create empty reactive list reacting to genelist input
rv <- reactiveValues(data = NULL, gene_lists = list())

# if genelist input is observed geneset tab is updated
observe({
  req(input$genelist)
  gene_data <- read_gene_list(input$genelist)
  rv$data <- gene_data
  rv$gene_lists[[input$genelist$name]] <- gene_data
  updateSelectInput(session, "genes",  choices = c("all genes", "kidney marker genes",
                                                   "tubular marker genes", names(rv$gene_lists)))
 
})

#--------------------------------------------#
# Reset
#-------------------------------------------#

 observeEvent(input$reset_test, {
   #sample(NULL)
   reset(id = "test", asis = FALSE)
   })

# observeEvent(input$reset, {
#   rv$data <- NULL
#   reset("genelist")
# })

#-----------------------------------------#
# select input: celltypes
#----------------------------------------#

# remove Gene.names column from displayed celltypes

names_ref <- reactive({
  colnames(reference())[colnames(reference()) != "Gene.names"]
})

names_test <- reactive({
  colnames(sample())[colnames(sample()) != "Gene.names"]
})

# add colnames of reference and test sample into select input box

observeEvent(reference(), {
  shinyWidgets::updatePickerInput(session, "reference_celltype_select",
                                  choices = names_ref(),
                                  selected = names_ref())
})

observeEvent(sample(), {
  shinyWidgets::updatePickerInput(session, "test_celltype_select",
                                  choices = names_test(),
                                  selected = names_test())
})

# define object containing selected celltypes and Gene.names

ref_to_display <- reactive(input[["reference_celltype_select"]])
test_to_display <- reactive(input[["test_celltype_select"]])

reference_select <- reactive({
  reference()[, c("Gene.names", input[["reference_celltype_select"]])]
})

test_select <- reactive({
  if(is.null(sample()[, c("Gene.names", input[["test_celltype_select"]])])) {
    return(NULL)
  } else {
     sample()[, c("Gene.names", input[["test_celltype_select"]])]
  }
})

#-------------------------------------------#
# selection of counts or cpm
#------------------------------------------#

selected_option <- reactive({
  input$counts
})

#------------------------------------------#
# reduce data, cor_data, ED_data, generate plots
#------------------------------------------#

# reduce data with selected celltypes
reduced_data <- reactive({
  if(is.null(rv_sample$data) & input$demo == "--") {
    return(NULL)
  } else if(is.null(rv_sample$data) & input$demo != "--") {
    red_data <- reduce_data(reference_select(), test_select(), geneset())
    return(red_data)
  } else if(!is.null(rv_sample$data) & input$demo == "--" & selected_option() == "counts") {
    red_data <- reduce_data_counts(reference_select(), test_select(), geneset())
    return(red_data)
  } else if(!is.null(rv_sample$data) & input$demo == "--" & selected_option() == "CPM") {
    red_data <- reduce_data(reference_select(), test_select(), geneset())
    return(red_data)
  }
})

# creates results table of cor and ED
results_cor <- reactive({
  if(is.null(rv_sample$data) & input$demo == "--") {
    return(NULL)
  } else {
    rho_data <- cor_data(reduced_data())
    return(rho_data)
  }
})

results_ED <- reactive({
  if(is.null(rv_sample$data) & input$demo == "--") {
  return(NULL)
} else {
  ed_data <- ED_data(reduced_data())
  return(ed_data)
}
})

# creates plot of cor and ED

output$cor_plot <- renderPlot({
  if(input$button > 0) {
    graph(results_cor())
  } else if(is.null(rv_sample$data) & input$demo == "--") {
    stop("Please select or upload test samples.")
  }
})

output$ED_plot <- renderPlot({
  if(input$button > 0) {
    graph_ED(results_ED())
  } else if(is.null(rv_sample$data) & input$demo == "--") {
    stop("Please select or upload test samples.")
  }
})


