#-----------------------------------------------#
#  Heatmap Code - ComplexHeatmap Version
#----------------------------------------------#
# Source the ComplexHeatmap function
source("functions/heatmap_complexheatmap.R")

# Load default marker genes
# Curated 40 tubular marker genes - most essential & well-known
tubular_markers_curated_40 <- c(
  # Proximal Tubule (PT) - 8 genes
  "SLC34A1",   # PT marker, phosphate transporter
  "LRP2",      # Megalin, classic PT marker
  "CUBN",      # Cubilin, PT endocytosis
  "SLC5A2",    # SGLT2, glucose transporter PTS1
  "SLC22A6",   # OAT1, organic anion transporter PTS1/2
  "SLC22A8",   # OAT3, PTS1
  "SLC13A3",   # Dicarboxylate transporter, PTS1/2/3
  "KAP",       # Kidney androgen-regulated protein, PTS2
  
  # Loop of Henle - Thin Limbs (DTL/ATL) - 6 genes
  "AQP1",      # Water channel, DTL
  "SLC14A2",   # Urea transporter, DTL
  "CLDN1",     # Claudin-1, tight junction, ATL
  "CLDN10",    # Claudin-10, ATL
  "CLCNKA",    # Chloride channel, ATL
  "PROX1",     # Transcription factor, ATL
  
  # Thick Ascending Limb (TAL) - 7 genes
  "UMOD",      # Uromodulin/Tamm-Horsfall, THE classic TAL marker
  "SLC12A1",   # NKCC2, sodium-potassium-chloride cotransporter
  "KCNJ1",     # ROMK, potassium channel
  "CASR",      # Calcium-sensing receptor
  "CLDN16",    # Claudin-16, CTAL
  "EGF",       # Epidermal growth factor, TAL/CTAL
  "TMEM207",   # Transmembrane protein, TAL
  
  # Distal Convoluted Tubule (DCT) - 7 genes
  "SLC12A3",   # NCC, THE classic DCT marker
  "TRPM6",     # Magnesium channel, DCT
  "CALB1",     # Calbindin D28k, DCT/CNT
  "PVALB",     # Parvalbumin, DCT1
  "TRPV5",     # Calcium channel, DCT2
  "WNK1",      # WNK kinase, DCT
  "WNK4",      # WNK4 kinase, DCT
  
  # Connecting Tubule (CNT) - 3 genes
  "ATP2B1",    # Calcium ATPase, CNT
  "KLK1",      # Kallikrein, CNT
  "SCNN1B",    # ENaC beta subunit, CNT/PC
  
  # Collecting Duct - Principal Cells (PC) - 5 genes
  "AQP2",      # Aquaporin 2, THE classic PC marker
  "AQP3",      # Aquaporin 3, PC
  "SCNN1G",    # ENaC gamma, sodium channel, PC
  "FXYD4",     # FXYD4, PC marker
  "AVPR2",     # Vasopressin V2 receptor, PC
  
  # Collecting Duct - Intercalated Cells (IC) - 4 genes
  "SLC4A1",    # AE1, bicarbonate exchanger, IC-A
  "ATP6V1B1",  # V-ATPase B1, proton pump, IC
  "SLC26A4",   # Pendrin, IC-B
  "FOXI1"      # Forkhead transcription factor, IC
)

# Use this in your heatmap code
markers <- data.frame(Gene.names = tubular_markers_curated_40)

# cell type assignment
gene_celltype <- read.csv("data/mg_all_celltype_assignment.csv")

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

# genesToPlot returns user selection OR NULL (which triggers default marker genes in heatmap function)
genesToPlot <- reactive({
  if(length(selected_genes$gene_list) > 0) {
    return(selected_genes$gene_list)
  } else {
    return(NULL)  # NULL means use default tubular marker genes
  }
})

# create reactive object containing Spearman results table ordered by descending rho
cor_order <- reactive({
  ordered_cor <- results_cor() %>%
    mutate(celltypes = ifelse(datatype == "Median of sample(s) vs. reference", paste0("ref_", celltypes),
                              ifelse(datatype == "Median of sample(s) vs. sample(s)", paste0("test_", celltypes), celltypes))) %>%
    arrange(desc(rho))
  return(ordered_cor)
})

# create ComplexHeatmap output
# In results_map reactive - add gene_celltype argument:
results_map <- reactive({
  req(input$button > 0)
  
  scale_bool <- input$scale == "scale by row"
  
  results_HM <- create_complex_heatmap_shiny(
    reference()[, c("Gene.names", input[["reference_celltype_select"]])],
    sample()[, c("Gene.names", input[["test_celltype_select"]])],
    markers,
    cor_order(),
    genesToPlot(),
    results_cor(),
    selected_option(),
    scale = scale_bool,
    gene_celltype = gene_celltype  # ADD THIS
  )
  return(results_HM)
})

# Reactive height based on number of genes
heatmap_height <- reactive({
  req(input$button > 0)
  n <- results_map()$n_genes
  calculate_heatmap_height(n)
})

# Render ComplexHeatmap with dynamic height
output$heatmap <- renderPlot({
  if(input$button > 0) {
    draw(results_map()$heatmap)
  }
}, height = function() heatmap_height())

