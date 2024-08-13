#----------------------------------------#
# load libraries
#---------------------------------------#

library(shiny)    
library(DT)
library(shinydashboard)
library(ggplot2)
library(readxl)
library(shinyWidgets)
library(shinyjs)

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
                               p("CellMatchR will run rank-based Spearman's correlations and Euclidean distance between your sample cells and reference cell types and
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
                                       choices = c("--","kidney primary cells - select cell types below", "HK-2 proximal tubule cell line", "mIMCD-3 cell line"))),
                  ),
                  fluidRow(
                    column(12,
                           selectInput("ref", 
                                       label = tags$label("2. Select a reference dataset",
                                                          actionButton("help_ref", label = NULL,icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;")),
                                       choices =  c("Ransick et al. (mouse, recommended)", "Park et al. (mouse)", "Lake et al. (human)", "Zhang et al. (human)" 
                                                    #,"Novella et al. (integration of multiple mouse references)"
                                                    )))
                  ),
                  br(),
                  br(),
                  fluidRow(
                    column(6,
                           selectInput("genes", 
                                       label = tags$label("3. Select a geneset...", 
                                                          actionButton("help_gene", label = NULL, icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;")),
                                                choices = c("all genes", "marker genes tubular celltypes"
                                                            #"marker genes all celltypes"
                                                            ))
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
                            30 cell clusters were identified containing epithelial, endothelial, stromal and immune cell types.</p>
                            <p>Data was downloaded from <a>https://github.com/qinzhu/kidneycellexplorer/tree/master/data</a> (2023/02/03). 
                            Detailed description of the cell types displayed in the legend was downloaded from: <a>https://cello.shinyapps.io/kidneycellexplorer/</a> (2023/02/20)</p>
                            <p>Reference: <em> Ransick, A., Lindström, N.O., Liu, J., Zhu, Q., Guo, J.-J., Alvarado, G.F., Kim, A.D., Black, H.G., Kim, J., McMahon, A.P., 2019. 
                            Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. Dev. Cell 51, 399-413.e7. <a>https://doi.org/10.1016/j.devcel.2019.10.005</a>
                            </em></p>" ),
                           hr(),
                           h4("Park et al., 2018"),
                           HTML("<p>First single cell RNA sequencing atlas of the mouse kidney generated with droplet-based single cell RNA sequencing. Kidneys were derived from 7 healthy male mice. 
                            24 cell clusters were identified containing epithelial, endothelial, stromal and immune cell types.</p>
                            <p>Data was downloaded from Gene Expression omnibus (GEO; accession no. GSE107585) and further processed with Seurat v2.3.4 for normalization. Cell types <q>novel 1</q> and <q>novel 2</q> were removed.</p>
                            <p>Reference: <em> Park, J., Shrestha, R., Qiu, C., Kondo, A., Huang, S., Werth, M., Li, M., Barasch, J., Suszták, K., 2018. 
                            Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease. Science 360, 758–763. <a>https://doi.org/10.1126/science.aar2131</a></em></p>"),
                           hr(),
                           h4("Kidney Precision Medicine Project (KPMP)"),
                           HTML("<p>Droplet-based single nucleus RNA sequencing dataset from human kidney generated with Chromium v3 platform. 58 kidney biopsies were derived from 35 healthy human donors.
                                  77 cell cluster were identified, including epithelial, endothelial, stromal, immune and neural cell types. </p>
                                  <p>The h5Seurat file of snRNA-seq data was downloaded from Kidney Cell Atlas website kpmp.org repository section on 2023/11/30. Cell type names and abbreviations were adapted from supplementary table 4 of <em>Lake et al.</em>.</p>
                                  <p>References: <br><em>https://www.kpmp.org accessed on 2023/11/30</em><br>
                                  <p><em>Lake et al., 2023. An atlas of healthy and injured cell states and niches 
                                  in the human kidney. Nature 619, 585–594.</em></p>")
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
                       h3("CellMatchR: Comparison of matching strategies", style = "color: #0072B2;"),
                       br(),
                       fluidRow(
                         column(6,
                          HTML("<p>CellMatchR statistically compares gene expression profiles derived from RNA sequencing data of whole tissues, primary cells and cell lines of the kidney to single cell and single nucleus RNA-seq references.
                          We tested multiple references, genesets and statistical approaches to evaluate which method provides the most accurate matching results. For this, we used bulk RNA-seq data from murine and human kidney primary cells
                          as positive controls and matched these to our references. We then assessed how often each method matches our positive control cell type to the correct corresponding reference cell type (figure 1).
                          Spearman's correlation, which compares the rank order of genes in between samples, performed slightly better than Euclidean distance which is based on the sum of pairwise gene-count (CPM) differences between samples. 
                          Reference Ransick <em>et al.</em> outperformed all other references.
                          Spearman's correlation using tubular marker gene expression of Ransick <em>et al.</em> against all positive controls matched 89% of the positive controls to the correct cell type.
                          Park <em>et al.</em> performed best using global gene expression and Spearman's correlation. Human snRNA-seq reference by Kidney Precision Medicine Project (KPMP) performed best using Euclidean distance and tubular marker gene expression. 16 out of 18
                          positive controls were derived from mice. KPMP might perform better as a reference dataset with human samples.</p>")
                          ),
                         column(6,
                                tags$div(
                                  tags$img(src = "score_graph.png", height = "100%", style = "max-width: 100%;"),
                                  tags$p(HTML("<b>Figure 1</b>&nbsp;<q>Positive controls</q> (n = 18) from bulk RNA-seq kidney primary cells and primary tissue were tested across 3 references, 2 algorithms and 3 gene sets. 
                                            Score of 100% means that CellMatchR always matched the correct reference cell type to the positive control cell type ."), 
                                       style = "font-size: 12px;"))
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
                           br()),
                        fluidRow(
                          column(6,
                                 "Service Project S1",
                                 br(),
                                 "Collaborative Research Center 1453",
                                 br(),
                                 "Nephrogenetics (NephGen)",
                                 br(),
                                 HTML("Website: <a>https://www.sfb1453.uni-freiburg.de/</a>")),
                          column(6,
                                 "Institute for Genetic Epidemiology",
                                 br(),
                                 "Universitätsklinikum Freiburg",
                                 br(),
                                 "Hugstetter Straße 49",
                                 br(),
                                 "79106 Freiburg, Germany",
                                 br(),
                                 HTML("Website: <a>https://www.uniklinik-freiburg.de/genetische-epidemiologie.html</a>"))
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
  
  #increase maximum upload size from 5MB to 10MB
  options(shiny.maxRequestSize = 20*1024^2)
  
  # source functions 
  source("functions/reduce_data.R")
  source("functions/cor_data.R")
  source("functions/graph_fct.R")
  source("functions/Euclidean_distance.R")
  source("functions/heatmap_function.R")
         
  
  # source dataset code
  source("shiny/server_all_datasets.R", local = TRUE)
  
  # source heatmap code
  source("shiny/heatmap_code.R", local = TRUE)
  
  # source help button code
  source("shiny/help_button_code.R", local = TRUE)
  
  # download code
  source("shiny/download_plots_tables.R", local = TRUE)
  
  # source reference box
  source("shiny/reference_box_server.R", local = TRUE)
  

  # create reactive legend object containing celltypes legend
  
  legend <- reactive({
    if(input$ref == "Ransick et al. (mouse, recommended)") {
      legend <- read.xlsx("data/Celltypes_Abbreviations.xlsx", sheet = "Ransick et al.")
    } else if(input$ref == "Park et al. (mouse)") {
      legend <- read.xlsx("data/Celltypes_Abbreviations.xlsx", sheet = "Park et al.")
    } else if(input$ref == "Zhang et al. (human)") {
      legend <- read.xlsx("data/Celltypes_Abbreviations.xlsx", sheet = "Zhang et al.")
    } else if(input$ref == "Lake et al. (human)") {
      legend <- read.xlsx("data/Celltypes_Abbreviations.xlsx", sheet = "Lake et al.")
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
        tabPanel("Euclidean distance",
                 fluidRow(actionLink("ED_interpret", "Click here to help with the results interpretation", style  = "padding: 25px;")),
                 #HTML("Results show <b>Euclidean distance</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b>.<br>
                               # Results are ordered from <b>lowest to highest Euclidean distance</b>. <b>Low Euclidean distance means higher similarity</b> between reference and sample gene expression.<br>
                              #<q>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Euclidean distance of all sample replicates.<br>")),
                 fluidRow(
                   column(8,
                          plotOutput("ED_plot", "600px", width = "100%", height = "800"),
                          downloadButton("DL_ED_plot", label = "Download plot", icon = icon("download"), style  = "padding: 0px;"),
                          downloadButton("DL_ED_table", label = "Download results as excel table", icon = icon("download"), style  = "padding: 0px;")),
                   # column(3,
                   #        div(style = "border: 1px solid #000; padding: 10px;",
                   #            HTML("Results show <b>Euclidean distance</b> of selected <b>reference cell type(s)</b> against the <b>median of all selected samples</b>.<br>
                   #            <q>Errorbars</q> indicate the <b>minimum</b> and <b>maximum</b> Euclidean distance of the samples."))))
        )),
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
                                                            actionButton("help_heatmap", label = NULL,  icon = icon("question-circle"), style = "font-size: 11px; border: none; background-color: transparent;")),
                                         choices =  c("no scaling", "scale by row")))
                   ),
                 br(),
                 fluidRow(
                   column(12, 
                          "By default, heatmap displays",
                          strong("kidney tubular marker genes"),
                          "of Lake et al. However, all genes of the datatable can be selected.",
                          p("Depending on the number of selected genes, not all gene names can be displayed on y-axis. Zoom-in for more details."),
                          hr())
                 ),
                 fluidRow(
                   column(12, plotlyOutput("heatmap", width = "100%", height = "800px")),
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



