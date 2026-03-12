#----------------------------------------#
# load libraries
#---------------------------------------#
library(shiny)    
library(DT)
library(shinydashboard)
library(ggplot2)
library("ggfortify")
library("ggrepel")
library("ggpubr")
library("ggdendro")
library("gridExtra")
library(readxl)
library(shinyWidgets)
library(shinyjs)
library(ComplexHeatmap)  
library(circlize)      
library(viridis)   
library(dplyr)
library(openxlsx)
library(writexl)

#------------------------------------------#
#  UI
#------------------------------------------#

ui <- dashboardPage( 
  title="NephGen | CellMatchR",
  skin = "black",
  dashboardHeader(
    title = span(
      tags$img(src = 'logo.png', height = '30'),
      )
  ),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "Home"),
      menuItem("Matching", tabName = "Matching"),
      menuItem("Reference Datasets", tabName = "References"),
      menuItem("About", tabName = "About")
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "Home", 
              column(12,
                     div(style = "background-color: white; padding: 20px;",
                         span(
                           style = "display: block; text-align: center;",
                           tags$img(src = 'logo_about.png', height = '80')),
                         hr(),
                         h3("How similar is your kidney cell line to actual kidney cell types?", style = "color: #0072B2;"),
                         #br(),
                         fluidRow(
                         column(9,
                                br(),
                                p("CellMatchR compares the genetic expression profile obtained from bulk RNA-sequencing of cell lines, primary cells or whole tissue to kidney single cell 
                                RNA-sequencing references."),
                                HTML("How to use CellMatchR:<il>
                                                    <li><b>Required</b>: bulk RNA-seq data containing raw counts or normalized to counts per million(CPM)
                                                    <li><b>Step 1</b>: Preprocess your bulk RNA-seq data to only include columns with gene identifiers and sample data</li>
                                                    <li><b>Step 2</b>: Upload formatted counts data</li>
                                                    <li><b>Step 3</b>: Select the reference publication, gene set and cell types you want to use</li>
                                                    </ul>"),
                                br(),
                               p("CellMatchR will run rank-based Spearman's correlations between your sample cells and reference cell types and
                              display the most similar cell type on top."),
                              p("Additionally, CellMatchR displays a heatmap to showcase gene expression profiles of selected cell types.")),
                         column(3,
                                tags$img(src = "nephron_culture_icon.png", height = "300"),
                                br(),
                                tags$p("created with Biorender.com", style = "font-size: 12px;")))
                         # fluidRow(
                         #   column(10, align = "center",
                         #          h3("Get Started"))
                         # )
              ))),
      tabItem(tabName = "Matching", 
              fluidRow(
                box(
                  title = "Input", 
                  status = "primary", solidHeader = TRUE, #height = 650,
                  #style = "max-height: 650px; overflow-y: auto;",
                  fluidRow(
                    useShinyjs(),
                    column(5, 
                           fileInput("test", 
                                     label = tags$label("1. Upload your sample...", 
                                                        actionButton("help_test", label = NULL,  icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;"),
                                                        actionLink("reset_test", "Reset", style  = "padding: 0px;")), 
                                     accept = c(".csv", ".xlsx"))
                           ),
                    column(1,
                           div(
                             style = "padding-top: 30px;margin-left: -20px;",
                             radioButtons("counts", label = NULL, c("counts", "CPM"), width = "100%") 
                           )),
                    column(6,
                           selectInput("demo", 
                                       label = tags$label("... or try out CellMatchR with our datasets:", 
                                                          actionButton("help_demo", label = NULL,icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;")),
                                       choices = c("--","kidney primary cells - select cell types below", "HK-2 proximal tubule cell line"))),
                  ),
                  fluidRow(
                    column(12,
                           selectInput("ref", 
                                       label = tags$label("2. Select a reference dataset",
                                                          actionButton("help_ref", label = NULL,icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;")),
                                       choices =  c("Ransick et al. (mouse, recommended)", "Park et al. (mouse)", "Lake et al. (human)"
                                                    , "Zhang et al. (human)" 
                                                    )))
                  ),
                  br(),
                  br(),
                  fluidRow(
                    column(6,
                           selectInput("genes", 
                                       label = tags$label("3. Select a geneset...", 
                                                          actionButton("help_gene", label = NULL, icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;")),
                                                choices = c("all genes", "kidney marker genes", "tubular marker genes")
                                       )
                                       ),
                    column(6,
                           fileInput("genelist", 
                                     label = tags$label(strong("... or upload your own marker gene list"), 
                                                        actionButton("help_upload", label = NULL,icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;"))))
                  ),
                  fluidRow(
                    column(5,
                           shinyWidgets::pickerInput(
                             "reference_celltype_select",
                             label = "4. Select celltypes from the reference",choices = c(""),
                             selected = c(""),
                             options = list("actions-box" = TRUE),
                             multiple = TRUE)),
                    column(1,
                           actionButton("help_select", label = NULL,  icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;")),
                    column(6,
                           shinyWidgets::pickerInput(
                             "test_celltype_select",
                             label = "Select celltypes from your sample(s)",
                             #label = tags$label("Select celltypes from your sample(s)",
                                                #actionButton("help_select", label = NULL,  icon = icon("question-circle"), style = "position: absolute; top: 0;font-size: 11px; border: none; background-color: transparent;padding: 0;margin: auto; display: inline; border-style: none;")),
                             choices = c(""),
                             selected = c(""),
                             options = list("actions-box" = TRUE),
                             multiple = TRUE)),
                  ),
                  br(),
                  fluidRow(
                    column(12, align = "center", actionButton("button", strong("Match!"), style = "background-color: #bd4861; border-color: black; color: white;", class = "btn-warning btn-lg"))),
                ),
                box(title = "Legend of selected reference",
                    status = "primary", solidHeader = TRUE, #height = 550,
                    collapsible = TRUE,
                    style = "height: 509px; overflow-y: auto;",  # adjusts height of box to previous tabbox, ensures scrolling
                    tableOutput("legend"))
                ),
                fluidRow(uiOutput("tabBox"))
                ),
      tabItem(tabName = "References",
              fluidRow(
                column(12,
                       div(style = "background-color: white; padding: 20px;",
                           h3("Which references were used?", style = "color: #0072B2;"),
                           hr(),
                           h4("Ransick et al., 2019"),
                           HTML("<p>Murine kidney cell atlas generated with Illumina Hi-Seq sequencing by 10X Genomics Chromium platform.
                            Kidneys were derived from 2 adult male and 2 adult female C57BL6/J mice. 
                            High resolution atlas as prior to cell dissociation kidney was subdivided into cortex, the outer medulla and the inner medulla. 
                            40 cell clusters were identified containing epithelial, endothelial, stromal and immune cell types.</p>
                            <p>Data was downloaded from <a>https://github.com/qinzhu/kidneycellexplorer/tree/master/data</a> (2023/02/03). 
                            Detailed description of the cell types displayed in the legend was downloaded from: <a>https://cello.shinyapps.io/kidneycellexplorer/</a> (2023/02/20)</p>
                            <p>Reference: <em> Ransick, A., Lindström, N.O., Liu, J., Zhu, Q., Guo, J.-J., Alvarado, G.F., Kim, A.D., Black, H.G., Kim, J., McMahon, A.P., 2019. 
                            Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. Dev. Cell 51, 399-413.e7. <a>https://doi.org/10.1016/j.devcel.2019.10.005</a>
                            </em></p>" ),
                           hr(),
                           h4("Park et al., 2018"),
                           HTML("<p>First single cell RNA sequencing atlas of the mouse kidney generated with droplet-based single cell RNA sequencing. Kidneys were derived from 7 healthy male mice. 
                            24 cell clusters were identified containing epithelial, endothelial, stromal and immune cell types. Cell type <q>Endo</q> was excluded.</p>
                            <p>Data was downloaded from Gene Expression omnibus (GEO; accession no. GSE107585) and further processed with Seurat v2.3.4 for normalization. Cell types <q>novel 1</q> and <q>novel 2</q> were removed.</p>
                            <p>Reference: <em> Park, J., Shrestha, R., Qiu, C., Kondo, A., Huang, S., Werth, M., Li, M., Barasch, J., Suszták, K., 2018. 
                            Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease. Science 360, 758–763. <a>https://doi.org/10.1126/science.aar2131</a></em></p>"),
                           hr(),
                           h4("Kidney Precision Medicine Project (KPMP)"),
                           HTML("<p>Droplet-based single-cell RNA sequencing based on 10x Genomics Chromium platform with Illumina Hi-Seq sequencing. 28 kidney biopsies were derived from 26 healthy human donors.
                                  77 cell cluster were identified, including epithelial, endothelial, stromal, immune and neural cell types. Only canonical cell types were used for analysis. </p>
                                  <p>The h5Seurat file was downloaded from Kidney Cell Atlas website kpmp.org repository section on 2025/04/08. Cell type names and abbreviations were adapted from supplementary table 4 of <em>Lake et al.</em>.</p>
                                  <p>References: <br><em>https://www.kpmp.org</em><br>
                                  <p><em>Lake et al., 2023. An atlas of healthy and injured cell states and niches 
                                  in the human kidney. Nature 619, 585–594.</em></p>"),
                           hr(),
                           h4("Zhang et al., 2021"),
                           HTML("<p>Droplet-based single cell RNA sequencing dataset of benign kidney and RCC tumor samples.
                                  26 cell cluster were identified, including epithelial, endothelial, stromal and immune cell types. 
                                  Only normal, non-tumor related cell types were used for the analysis. Cell types <q>ua</q>, <q>UC</q> and <q>unknown</q> were excluded.</p>
                                  <p>Cell type assignments and RNA-seq count matrix were downloaded from 
                                  Gene Expression Omnibus (GEO; accession no. GSE159115, 2024/10/03).</p>
                                  <p>Reference: <em>Zhang et al., 2021. Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response. 
                                   Proc Natl Acad Sci U S A.</em></p>")
                       )
                       ))),
      tabItem(tabName = "About",
              fluidRow(
                    column(10,
                       div(style = "background-color: white; padding: 20px;",
                       span(
                         style = "display: block; text-align: center;",
                         tags$img(src = 'logo_about.png', height = '80')),
                       hr(),
                       h3("CellMatchR Publication", style = "color: #0072B2;"),
                       br(),
                       fluidRow(
                         column(6,
                                HTML('<p>For more information on the development and usage of this tool, refer to our <a href="[DOI_URL]" target="_blank">publication</a> (under review).
                                Detailed documentation of the code is available on <a href="https://github.com/genepi-freiburg/CellMatchR/tree/main" target="_blank">GitHub</a>.
                                Comparison of bulk RNA-seq data to scRNA-seq kidney references can also be performed using the TabPFN machine learning model — detailed instructions are available 
                                <a href="https://github.com/genepi-freiburg/CellMatchR/tree/main/TabPFN_scripts" target="_blank">here</a>.</p>')
                          )
                       ))
                    )),
              fluidRow(
                column(10,
                       title = "Legal notice(Impressum)",
                       div(
                        style = "background-color: white; padding: 20px;",
                        h3("Legal notice (Impressum)", style = "color: #0072B2;"),
                        hr(),
                        h4("Responsible for the content of this website:",
                           br(),
                           br()),
                        fluidRow(
                          column(6,
                                 img(src = "NephGen_Logo.svg", height = "100px", style = "margin-bottom: 10px;"),
                                 br(),
                                 "Service Project S1",
                                 br(),
                                 "Collaborative Research Center 1453",
                                 br(),
                                 "Nephrogenetics (NephGen)",
                                 br(),
                                 HTML("Website: <a href='https://www.sfb1453.uni-freiburg.de/' target='_blank'>https://www.sfb1453.uni-freiburg.de/</a>")
                                 ),
                          column(6,
                                 img(src = "EPI_logo_full.svg", height = "100px", style = "margin-bottom: 10px;"),
                                 br(),
                                 "Institute of Epidemiology and Prevention",
                                 br(),
                                 "Universitätsklinikum Freiburg",
                                 br(),
                                 "Hugstetter Straße 49",
                                 br(),
                                 "79106 Freiburg, Germany",
                                 br(),
                                 HTML("Website: <a href='https://www.uniklinik-freiburg.de/epidemiologie.html' target='_blank'>https://www.uniklinik-freiburg.de/epidemiologie.html</a>")
                          )
                        ),
                        br(),
                        HTML("<u>We are happy to receive feedback. For this, please contact:</u>"),
                        br(),
                        fluidRow(
                          column(6,
                                 "Mona Schoberth: mona.schoberth@uniklinik-freiburg.de"),
                          column(6,
                                 "Dr. Stefan Haug:  stefan.haug@uniklinik-freiburg.de"))))),
              fluidRow(
                column(10,
                       title = "Legal notice(Impressum)",
                       div(
                         style = "background-color: white; padding: 20px;",
                         h3("Disclaimer", style = "color: #0072B2;"),
                         hr(),
                         p("The	author	assumes	no	responsibility	for	the	topicality,	correctness,	completeness	or	quality	of information	provided.	Liability	claims	against	the	author	which	relate	to	material	or	immaterial	
                          nature	caused	by	the	use	or	misuse	of	any	information	provided	through	the	use	of	incorrect	or incomplete	information	are	excluded	unless	the	author	is	not	intentional	or	grossly	negligent	fault.	
                          The	author	reserves	the	right	to	change	parts	of	the	site	without	prior	notice,	add	to,	delete	or	cease	publication	temporarily	or	permanently."))
              )
      )
    )
  )
))



#-----------------------------------------------#
#  server code
#-----------------------------------------------#

server <- function(input,output, session) {
  
  #increase maximum upload size from 5MB to 100MB
  options(shiny.maxRequestSize = 100*1024^2)
  
  # source functions 
  source("functions/reduce_data.R")
  source("functions/cor_data.R")
  source("functions/graph_fct.R")
 # source("functions/Euclidean_distance.R")
  source("functions/heatmap_complexheatmap.R")
         
  
  # source dataset code
  source("shiny/server_all_datasets.R", local = TRUE)
  
  # source heatmap code
  source("shiny/complexheatmap_code.R", local = TRUE)
  
  # source help button code
  source("shiny/help_button_code.R", local = TRUE)
  
  # download code
  source("shiny/download_plots_tables.R", local = TRUE)
  
  # source reference box
  source("shiny/reference_box_server.R", local = TRUE)
  

  # create reactive legend object containing celltypes legend
  
  legend <- reactive({
    if(input$ref == "Ransick et al. (mouse, recommended)") {
      legend <- read.xlsx("data/Abbreviations_Celltypes.xlsx", sheet = "Ransick et al.")
    } else if(input$ref == "Park et al. (mouse)") {
      legend <- read.xlsx("data/Abbreviations_Celltypes.xlsx", sheet = "Park et al.")
    } else if(input$ref == "Zhang et al. (human)") {
      legend <- read.xlsx("data/Abbreviations_Celltypes.xlsx", sheet = "Zhang et al.")
    } else if(input$ref == "Lake et al. (human)") {
      legend <- read.xlsx("data/Abbreviations_Celltypes.xlsx", sheet = "Lake et al.")
    }
    legend
  })
  
   output$legend <- renderTable(
       legend()
   )
  
  # create reference box

  output$ReferenceBox <- renderUI({
    req(referenceContent())
    
    HTML(referenceContent())
    })
  
  # create results box upon pressing Match button
   
  output$tabBox <- renderUI({
    if(input$button == 0 || is.null(rv_sample$data) & input$demo == "--" ) return(NULL)
    
    fluidRow(
      ## hide all output errors in Shiny for the user
      
      tags$style(type="text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"
      ),
      #style = "width: 3000px; padding: 18px; overflow-y: auto;",
      div(
        style = "padding: 18px;",
        tabBox(
          title = "Results",
          width = 12,
          tabPanel("Spearman correlation", 
                 fluidRow(actionLink("rho_interpret", "Click here to help with the results interpretation", style  = "padding: 25px;")),# HTML("Results show <b>Spearman's rho</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b>.<br>
                            #Results are ordered from <b>highest to lowest rho</b>. <b>High rho means higher correlation</b> between reference and sample gene expression ranks.<br>
                            #<q>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Spearman's rho of sample replicates.<br>")),
                 fluidRow(
                   column(width = 8,
                          plotOutput("cor_plot", "600px", width = "100%", height = "800"),
                          downloadButton("DL_cor_plot", label = "Download plot", icon = icon("download"), style  = "padding: 0px;"),
                          downloadButton("DL_cor_table", label = "Download results as excel table", icon = icon("download"), style  = "padding: 0px; padding-left: 1cm;"))
                  # column(3,
                  #        div(style = "border: 1px solid #000; padding: 10px;",
                  #           HTML("Results show <b>Spearman's rho</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b>.<br>
                  #           <q>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Spearman's rho of the samples."))))
                 )),
        # tabPanel("Euclidean distance",
        #          fluidRow(actionLink("ED_interpret", "Click here to help with the results interpretation", style  = "padding: 25px;")),
        #          #HTML("Results show <b>Euclidean distance</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b>.<br>
        #                        # Results are ordered from <b>lowest to highest Euclidean distance</b>. <b>Low Euclidean distance means higher similarity</b> between reference and sample gene expression.<br>
        #                       #<q>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Euclidean distance of all sample replicates.<br>")),
        #          fluidRow(
        #            column(8,
        #                   plotOutput("ED_plot", "600px", width = "100%", height = "800"),
        #                   downloadButton("DL_ED_plot", label = "Download plot", icon = icon("download"), style  = "padding: 0px;"),
        #                   downloadButton("DL_ED_table", label = "Download results as excel table", icon = icon("download"), style  = "padding: 0px;")),
        #            # column(3,
        #            #        div(style = "border: 1px solid #000; padding: 10px;",
        #            #            HTML("Results show <b>Euclidean distance</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b>.<br>
        #            #            <q>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Euclidean distance of the samples."))))
        # )),
        tabPanel("Heatmap",
                 fluidRow(
                   column(width = 4,
                          selectizeInput("expression_genes_input",
                                         label = tags$div("Gene search", HTML("&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"),
                                                          actionLink("reset_genes", "Reset", style  = "padding: 0px;")
                                         ),
                                         choices = NULL,
                                         options = list(create = TRUE, placeholder = "Choose or enter gene name(s) here"),
                                         multiple= TRUE)
                   ),
                   column(4, selectInput("scale",
                                         label = tags$label("Heatmap options",
                                                            actionButton("help_heatmap", label = NULL,  icon = icon("question-circle"), 
                                                                         style = "font-size: 11px; border: none; background-color: transparent;")),
                                         choices =  c("no scaling", "scale by row")))
                   ),
                 br(),
                 fluidRow(
                   column(12, 
                          "By default, heatmap displays",
                          strong(" 40 selected tubular marker genes."),
                          "However, all genes of the datatable can be selected.",
                          p("Genes are labelled with the respective cell types and appended with a dot and an incremental number in case they are markers for multiple tubular cell types"),
                          hr())
                 ),
                 fluidRow(
                   column(12, plotOutput("heatmap", width = "100%", height = "auto")),
                 ),
                 fluidRow(
                   column(4,
                          downloadButton("DL3", label = "Download heatmap datatable", icon = icon("download"))))
        ))
       )
      #,
      # box(
      #   title = "Legend",
      #   style = "height: 800px; overflow-y: auto;",  # adjusts height of box to previous tabbox, ensures scrolling
      #   tableOutput("legend")
      # )
    )
    
  })
}

#-------------------------------------------------#
#  Start App
#------------------------------------------------#

shinyApp(ui, server)



