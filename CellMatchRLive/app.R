# CellMatchR - Shinylive Version
# Cell type matching using Spearman correlation
# Deployable to GitHub Pages via Shinylive

library(shiny)
library(shinydashboard)
library(ggplot2)
library(shinyWidgets)
library(shinyjs)
library(dplyr)
library(viridisLite)
library(tidyr)

# ============================================================================
# Curated 40 tubular marker genes
# ============================================================================

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

curated_markers_df <- data.frame(Gene.names = tubular_markers_curated_40)

# ============================================================================
# Core Functions
# ============================================================================

reduce_data <- function(ref, test, marker_genes = NULL) {
  test$Gene.names <- toupper(test$Gene.names)
  ref$Gene.names <- toupper(ref$Gene.names)

  ref <- ref[!duplicated(ref$Gene.names), ]
  test <- test[!duplicated(test$Gene.names), ]

  test_filt <- na.omit(test)
  ref_filt <- na.omit(ref)

  if (is.null(marker_genes)) {
    gene_overlap <- data.frame(Gene.names = intersect(test_filt$Gene.names, ref_filt$Gene.names))
  } else {
    gene_overlap <- intersect(test_filt$Gene.names, ref_filt$Gene.names)
    marker_genes <- toupper(marker_genes$Gene.names)
    gene_overlap <- data.frame(Gene.names = intersect(gene_overlap, marker_genes))
  }

  test_filt <- dplyr::left_join(gene_overlap, test, by = "Gene.names")
  ref_filt <- dplyr::left_join(gene_overlap, ref, by = "Gene.names")

  rownames(test_filt) <- test_filt$Gene.names
  test_filt$Gene.names <- NULL

  rownames(ref_filt) <- ref_filt$Gene.names
  ref_filt$Gene.names <- NULL

  test_median <- apply(test_filt, 1, median)

  QC <- round((nrow(test_filt) / nrow(test)) * 100, 2)
  QC <- paste0("(", QC, "% of ", nrow(test), " genes used)")

  list(
    ref_filt = ref_filt,
    test_filt = test_filt,
    test_median = as.numeric(test_median),
    QC = QC,
    n_genes = nrow(test_filt)
  )
}

reduce_data_counts <- function(ref, test, marker_genes = NULL) {
  test <- na.omit(test)

  cpm_test <- as.data.frame(apply(
    test[, -grep("Gene.names", colnames(test)), drop = FALSE],
    MARGIN = 2,
    FUN = function(x) x / sum(x) * 1e6
  ))
  cpm_test$Gene.names <- test$Gene.names
  test <- cpm_test

  reduce_data(ref, test, marker_genes)
}

cor_data <- function(reduced_data) {
  ref_test <- cbind(
    data.frame(reduced_data$ref_filt),
    data.frame(reduced_data$test_filt),
    test_median = reduced_data$test_median
  )
  colnames(ref_test)[(ncol(reduced_data$ref_filt) + 1):(ncol(reduced_data$ref_filt) + ncol(reduced_data$test_filt))] <-
    colnames(reduced_data$test_filt)

  rho <- as.data.frame(cor(ref_test, method = "spearman"))
  rho <- rho[-nrow(rho), ]

  n_ref <- ncol(reduced_data$ref_filt)
  n_test <- ncol(reduced_data$test_filt)

  min_rho <- apply(rho[, (n_ref + 1):(n_ref + n_test + 1)], 1, min)
  max_rho <- apply(rho[, (n_ref + 1):(n_ref + n_test + 1)], 1, max)

  results <- data.frame(
    rho = rho$test_median,
    celltypes = rownames(rho),
    datatype = ifelse(
      seq_len(nrow(rho)) <= n_ref,
      "Median of sample(s) vs. reference",
      "Median of sample(s) vs. sample(s)"
    ),
    min = min_rho,
    max = max_rho,
    row.names = NULL
  )

  dups <- duplicated(results$celltypes)
  results$celltypes[dups] <- paste0(" ", results$celltypes[dups])

  results
}

# ED_data <- function(reduced_data) {
#   # Euclidean distance - commented out for Shinylive version
# }

graph <- function(results) {
  results$min[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA
  results$max[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA

  ggplot(results, aes(x = rho, y = reorder(celltypes, rho), fill = datatype)) +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%.2f", rho)), vjust = 0.5, hjust = -0.2, size = 4.2, color = "black") +
    geom_errorbar(aes(xmin = min, xmax = max), width = 0.2, linewidth = 0.5, color = "grey40") +
    geom_vline(xintercept = 0.8, linetype = "dashed", colour = "black") +
    geom_vline(xintercept = 0.6, linetype = "dashed", colour = "grey60") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2), expand = expansion(mult = c(0, 0.05))) +
    expand_limits(x = c(0, 1.05 * max(results$rho, na.rm = TRUE))) +
    facet_grid(datatype ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Spearman's rho", y = NULL) +
    scale_fill_manual(values = c("#B85C00", "#3B75A3")) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text.y = element_blank()
    )
}

# graph_ED <- function(results) {
#   # Euclidean distance plot - commented out for Shinylive version
# }

# Heatmap function
normalize_expression_by_row <- function(x) {
  if (max(x) == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - min(x)) / (max(x) - min(x)))
  }
}

calculate_heatmap_height <- function(n_genes) {
  base_height <- 200
  height_per_gene <- 18
  min_height <- 400
  max_height <- 2000
  calculated <- base_height + (n_genes * height_per_gene)
  max(min_height, min(max_height, calculated))
}

heatmap_fct <- function(reduced_data, cor_results, markers = NULL, genesToPlot = NULL,
                        scale = FALSE, gene_celltype = NULL) {
  data_filt <- reduced_data

  # Filter genes: user-selected > curated markers > all marker genes (cap at 100)
  if (!is.null(genesToPlot) && length(genesToPlot) > 0) {
    data_filt$ref_filt <- data_filt$ref_filt[rownames(data_filt$ref_filt) %in% genesToPlot, , drop = FALSE]
    data_filt$test_filt <- data_filt$test_filt[rownames(data_filt$test_filt) %in% genesToPlot, , drop = FALSE]
  } else if (!is.null(markers)) {
    marker_gene_list <- markers$Gene.names
    if (length(marker_gene_list) > 100) {
      marker_gene_list <- marker_gene_list[1:100]
    }
    data_filt$ref_filt <- data_filt$ref_filt[rownames(data_filt$ref_filt) %in% marker_gene_list, , drop = FALSE]
    data_filt$test_filt <- data_filt$test_filt[rownames(data_filt$test_filt) %in% marker_gene_list, , drop = FALSE]
  }

  # Prefix columns
  ref_filt <- data_filt$ref_filt
  colnames(ref_filt) <- paste0("ref_", colnames(ref_filt))

  test_filt <- data_filt$test_filt
  colnames(test_filt) <- paste0("test_", colnames(test_filt))

  # Order columns by descending rho: test on left, ref on right
  cor_order <- cor_results %>%
    dplyr::mutate(celltypes_prefixed = ifelse(
      datatype == "Median of sample(s) vs. reference",
      paste0("ref_", trimws(celltypes)),
      paste0("test_", trimws(celltypes))
    )) %>%
    dplyr::arrange(dplyr::desc(rho))

  test_col_order <- cor_order$celltypes_prefixed[cor_order$celltypes_prefixed %in% colnames(test_filt)]
  ref_col_order <- cor_order$celltypes_prefixed[cor_order$celltypes_prefixed %in% colnames(ref_filt)]

  data <- cbind(
    test_filt[, test_col_order, drop = FALSE],
    ref_filt[, ref_col_order, drop = FALSE]
  )

  # Cell type grouping of genes
  use_facets <- !is.null(gene_celltype) && (is.null(genesToPlot) || length(genesToPlot) == 0)

  if (use_facets) {
    present_genes <- rownames(data)
    gene_ct_filt <- gene_celltype[gene_celltype$Gene.names %in% present_genes, ]

    expanded_rows <- list()
    gene_labels <- c()
    cell_type_labels <- c()

    for (gene in present_genes) {
      cts <- gene_ct_filt$abbreviation[gene_ct_filt$Gene.names == gene]
      if (length(cts) == 0) {
        expanded_rows[[length(expanded_rows) + 1]] <- data[gene, , drop = FALSE]
        gene_labels <- c(gene_labels, gene)
        cell_type_labels <- c(cell_type_labels, "Other")
      } else if (length(cts) == 1) {
        expanded_rows[[length(expanded_rows) + 1]] <- data[gene, , drop = FALSE]
        gene_labels <- c(gene_labels, gene)
        cell_type_labels <- c(cell_type_labels, cts)
      } else {
        for (j in seq_along(cts)) {
          expanded_rows[[length(expanded_rows) + 1]] <- data[gene, , drop = FALSE]
          gene_labels <- c(gene_labels, paste0(gene, ".", j))
          cell_type_labels <- c(cell_type_labels, cts[j])
        }
      }
    }

    data_expanded <- do.call(rbind, expanded_rows)
    rownames(data_expanded) <- gene_labels

    sort_order <- order(cell_type_labels)
    data_expanded <- data_expanded[sort_order, ]
    cell_type_labels <- cell_type_labels[sort_order]

    data <- data_expanded
  } else {
    cell_type_labels <- NULL
  }

  # Log2 transformation
  data_log2 <- log2(data + 1)

  if (scale) {
    scale_data <- as.data.frame(t(apply(data_log2, 1, normalize_expression_by_row)))
    colnames(scale_data) <- colnames(data_log2)
    data_log2 <- scale_data
  }

  # Column annotation: Sample/Reference grouping
  n_test <- length(test_col_order)
  n_ref <- length(ref_col_order)
  all_cols <- c(test_col_order, ref_col_order)

  # Dynamic font size based on gene count
  n_genes <- nrow(data_log2)
  row_fontsize <- if (n_genes > 100) 6 else if (n_genes > 50) 9 else 11

  # Row ordering: clustering for user-selected genes, cell-type order for default
  if (use_facets) {
    gene_order <- rownames(data_log2)
  } else if (n_genes > 1) {
    hc <- hclust(dist(data_log2), method = "ward.D2")
    gene_order <- rownames(data_log2)[hc$order]
  } else {
    gene_order <- rownames(data_log2)
  }

  # Reshape to long format for ggplot
  plot_df <- data_log2
  plot_df$gene <- rownames(plot_df)
  if (use_facets) {
    plot_df$cell_type <- cell_type_labels
  }
  plot_long <- tidyr::pivot_longer(plot_df,
    cols = all_of(all_cols),
    names_to = "column", values_to = "value"
  )
  plot_long$gene <- factor(plot_long$gene, levels = rev(gene_order))
  plot_long$column <- factor(plot_long$column, levels = all_cols)

  # Insert spacer column between Sample and Reference
  spacer_name <- "   "
  spacer_df <- expand.grid(
    gene = levels(plot_long$gene),
    column = spacer_name,
    stringsAsFactors = FALSE
  )
  spacer_df$value <- NA
  if (use_facets) {
    gene_ct_map <- setNames(cell_type_labels, rownames(data_log2))
    spacer_df$cell_type <- gene_ct_map[spacer_df$gene]
  }
  spacer_df$gene <- factor(spacer_df$gene, levels = levels(plot_long$gene))

  col_levels <- c(test_col_order, spacer_name, ref_col_order)
  plot_long$column <- factor(plot_long$column, levels = col_levels)
  spacer_df$column <- factor(spacer_df$column, levels = col_levels)
  plot_all <- rbind(plot_long, spacer_df)

  # Color column labels: orange for Sample, blue for Reference
  col_label_colors <- ifelse(
    col_levels %in% test_col_order, "#D55E00",
    ifelse(col_levels %in% ref_col_order, "#0072B2", "white")
  )

  # Build heatmap
  p <- ggplot(plot_all, aes(x = column, y = gene, fill = value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colors = viridisLite::plasma(100),
      na.value = "white",
      name = if (scale) "Scaled\nexpression" else "log2(x+1)"
    ) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0.5, size = 9,
                                 color = col_label_colors, face = "bold"),
      axis.text.y = element_text(size = row_fontsize, color = "black"),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      legend.key.height = unit(1.5, "cm"),
      plot.margin = margin(5, 5, 5, 5)
    )

  if (use_facets) {
    p <- p +
      facet_grid(cell_type ~ ., scales = "free_y", space = "free_y", switch = "y") +
      theme(
        strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.5),
        strip.text.y.left = element_text(angle = 0, size = 9, face = "bold", margin = margin(2, 4, 2, 4)),
        strip.placement = "outside",
        panel.spacing.y = unit(2, "pt")
      )
  }

  gene_names <- rownames(data_log2)
  download_matrix <- cbind(gene_names, data_log2)

  list(plot = p, download_matrix = download_matrix, n_genes = n_genes)
}

# ============================================================================
# Cell Type Legends - Read from CSV files (generated from Excel)
# ============================================================================

# Load legend CSV files
legends <- list(
  "Ransick et al. (mouse, recommended)" = read.csv("data/legend_Ransick_et_al.csv"),
  "Park et al. (mouse)" = read.csv("data/legend_Park_et_al.csv"),
  "Lake et al. (human)" = read.csv("data/legend_Lake_et_al.csv"),
  "Zhang et al. (human)" = read.csv("data/legend_Zhang_et_al.csv")
)

# ============================================================================
# UI - Using shinydashboard (same as original)
# ============================================================================

ui <- dashboardPage(
  title = "NephGen | CellMatchR",
  skin = "black",
  dashboardHeader(
    title = span(tags$img(src = "logo.png", height = "30"))
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
    tags$head(
      tags$style(HTML("
        .shinylive-note { background-color: #fff3cd; padding: 15px; border-radius: 5px; border-left: 4px solid #ffc107; margin: 15px 0; }
      "))
    ),
    tabItems(
      # Home Tab
      tabItem(
        tabName = "Home",
        column(12,
          div(style = "background-color: white; padding: 20px;",
            span(style = "display: block; text-align: center;",
                 tags$img(src = "logo_about.png", height = "80")),
            hr(),
            h3("How similar is your kidney cell line to actual kidney cell types?", style = "color: #0072B2;"),
            fluidRow(
              column(9,
                br(),
                p("CellMatchR compares the genetic expression profile obtained from bulk RNA-sequencing of cell lines, primary cells or whole tissue to kidney single cell RNA-sequencing references."),
                HTML("How to use CellMatchR:<ul>
                  <li><b>Required</b>: bulk RNA-seq data containing raw counts or normalized to counts per million (CPM)</li>
                  <li><b>Step 1</b>: Preprocess your bulk RNA-seq data to only include columns with gene identifiers and sample data</li>
                  <li><b>Step 2</b>: Upload formatted counts data (CSV format)</li>
                  <li><b>Step 3</b>: Select the reference publication, gene set and cell types you want to use</li>
                </ul>"),
                br(),
                p("CellMatchR will run rank-based Spearman's correlations between your sample cells and reference cell types and display the most similar cell type on top."),
                p("Additionally, CellMatchR displays a heatmap to showcase gene expression profiles of selected cell types."),
                br(),
                div(class = "shinylive-note",
                  HTML("<b>Note:</b> This Shinylive version runs entirely in your browser using WebAssembly.
                       Initial load may take 30-60 seconds as the R runtime is downloaded.
                       All data processing happens locally - <b>your data is never uploaded to a server</b>."))
              ),
              column(3,
                tags$img(src = "nephron_culture_icon.png", height = "300"),
                br(),
                tags$p("created with Biorender.com", style = "font-size: 12px;")
              )
            )
          )
        )
      ),

      # Matching Tab
      tabItem(
        tabName = "Matching",
        fluidRow(
          box(
            title = "Input",
            status = "primary",
            solidHeader = TRUE,
            fluidRow(
              useShinyjs(),
              column(5,
                fileInput("test",
                          label = tags$label("1. Upload your sample (CSV)",
                                             actionLink("reset_test", "Reset", style = "padding: 0px;")),
                          accept = ".csv")
              ),
              column(1,
                div(style = "padding-top: 30px; margin-left: -20px;",
                    radioButtons("counts", label = NULL, c("counts", "CPM"), width = "100%"))
              ),
              column(6,
                selectInput("demo",
                            label = "... or try out CellMatchR with our datasets:",
                            choices = c("--", "kidney primary cells - select cell types below",
                                        "HK-2 proximal tubule cell line", "mIMCD-3 cell line"))
              )
            ),
            fluidRow(
              column(12,
                selectInput("ref",
                            label = "2. Select a reference dataset",
                            choices = c("Ransick et al. (mouse, recommended)", "Park et al. (mouse)",
                                        "Lake et al. (human)", "Zhang et al. (human)"))
              )
            ),
            fluidRow(
              column(12,
                fileInput("custom_ref", "... or upload your own reference (CSV)", accept = ".csv")
              )
            ),
            br(),
            fluidRow(
              column(6,
                selectInput("genes",
                            label = "3. Select a geneset...",
                            choices = c("all genes", "kidney marker genes", "tubular marker genes"))
              ),
              column(6,
                fileInput("genelist",
                          label = tags$label(strong("... or upload your own marker gene list")),
                          accept = ".csv")
              )
            ),
            fluidRow(
              column(6,
                pickerInput("reference_celltype_select",
                            label = "4. Select celltypes from the reference",
                            choices = c(""),
                            selected = c(""),
                            options = list(`actions-box` = TRUE),
                            multiple = TRUE)
              ),
              column(6,
                pickerInput("test_celltype_select",
                            label = "Select celltypes from your sample(s)",
                            choices = c(""),
                            selected = c(""),
                            options = list(`actions-box` = TRUE),
                            multiple = TRUE)
              )
            ),
            br(),
            fluidRow(
              column(12, align = "center",
                     actionButton("button", strong("Match!"),
                                  style = "background-color: #bd4861; border-color: black; color: white;",
                                  class = "btn-warning btn-lg"))
            )
          ),
          box(
            title = "Legend of selected reference",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            style = "height: 509px; overflow-y: auto;",
            tableOutput("legend")
          )
        ),
        fluidRow(uiOutput("tabBox"))
      ),

      # References Tab
      tabItem(
        tabName = "References",
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
                   <p>Data was downloaded from <a href='https://github.com/qinzhu/kidneycellexplorer/tree/master/data' target='_blank'>github.com/qinzhu/kidneycellexplorer</a> (2023/02/03).
                   Detailed description of the cell types displayed in the legend was downloaded from: <a href='https://cello.shinyapps.io/kidneycellexplorer/' target='_blank'>cello.shinyapps.io/kidneycellexplorer</a> (2023/02/20)</p>
                   <p>Reference: <em>Ransick, A., Lindstr&ouml;m, N.O., Liu, J., Zhu, Q., Guo, J.-J., Alvarado, G.F., Kim, A.D., Black, H.G., Kim, J., McMahon, A.P., 2019.
                   Single-Cell Profiling Reveals Sex, Lineage, and Regional Diversity in the Mouse Kidney. Dev. Cell 51, 399-413.e7.
                   <a href='https://doi.org/10.1016/j.devcel.2019.10.005' target='_blank'>https://doi.org/10.1016/j.devcel.2019.10.005</a></em></p>"),
              hr(),
              h4("Park et al., 2018"),
              HTML("<p>First single cell RNA sequencing atlas of the mouse kidney generated with droplet-based single cell RNA sequencing.
                   Kidneys were derived from 7 healthy male mice. 24 cell clusters were identified containing epithelial, endothelial, stromal and immune cell types.
                   Cell type &lsquo;Endo&rsquo; was excluded.</p>
                   <p>Data was downloaded from Gene Expression omnibus (GEO; accession no. GSE107585) and further processed with Seurat v2.3.4 for normalization.
                   Cell types &lsquo;novel 1&rsquo; and &lsquo;novel 2&rsquo; were removed.</p>
                   <p>Reference: <em>Park, J., Shrestha, R., Qiu, C., Kondo, A., Huang, S., Werth, M., Li, M., Barasch, J., Suszt&aacute;k, K., 2018.
                   Single-cell transcriptomics of the mouse kidney reveals potential cellular targets of kidney disease. Science 360, 758-763.
                   <a href='https://doi.org/10.1126/science.aar2131' target='_blank'>https://doi.org/10.1126/science.aar2131</a></em></p>"),
              hr(),
              h4("Kidney Precision Medicine Project (KPMP)"),
              HTML("<p>Droplet-based single cell RNA sequencing dataset from human kidney generated with Chromium v3 platform.
                   Only the subset of cells from healthy donors (n = 26) was used.
                   77 cell clusters were identified, including epithelial, endothelial, stromal, immune and neural cell types.
                   Only canonical cell types were used for analysis.</p>
                   <p>The h5Seurat file was downloaded from Kidney Cell Atlas website <a href='https://www.kpmp.org' target='_blank'>kpmp.org</a> repository section on 2025/04/08.
                   Cell type names and abbreviations were adapted from supplementary table 4 of <em>Lake et al.</em></p>
                   <p>References:<br><em><a href='https://www.kpmp.org' target='_blank'>https://www.kpmp.org</a> accessed on 2023/11/30</em></p>
                   <p><em>Lake et al., 2023. An atlas of healthy and injured cell states and niches in the human kidney. Nature 619, 585-594.
                   <a href='https://doi.org/10.1038/s41586-023-05769-3' target='_blank'>https://doi.org/10.1038/s41586-023-05769-3</a></em></p>"),
              hr(),
              h4("Zhang et al., 2021"),
              HTML("<p>Droplet-based single cell RNA sequencing dataset of benign kidney and RCC tumor samples.
                   26 cell clusters were identified, including epithelial, endothelial, stromal and immune cell types.
                   Only normal, non-tumor related cell types were used for the analysis.
                   Cell types &lsquo;ua&rsquo;, &lsquo;UC&rsquo; and &lsquo;unknown&rsquo; were excluded.</p>
                   <p>Cell type assignments and RNA-seq count matrix were downloaded from
                   Gene Expression Omnibus (GEO; accession no. GSE159115, 2024/10/03).</p>
                   <p>Reference: <em>Zhang et al., 2021. Single-cell analyses of renal cell cancers reveal insights into tumor microenvironment, cell of origin, and therapy response.
                   Proc Natl Acad Sci U S A.</em></p>")
            )
          )
        )
      ),

      # About Tab
      tabItem(
        tabName = "About",
        fluidRow(
          column(10,
            div(style = "background-color: white; padding: 20px;",
              span(style = "display: block; text-align: center;",
                   tags$img(src = "logo_about.png", height = "80")),
              hr(),
              h3("CellMatchR: Comparison of matching strategies", style = "color: #0072B2;"),
              br(),
              fluidRow(
                column(6,
                  HTML("<p>CellMatchR statistically compares gene expression profiles derived from RNA sequencing data
                       of whole tissues, primary cells and cell lines of the kidney to single cell and single nucleus RNA-seq references.
                       We tested multiple references, genesets and statistical approaches to evaluate which method provides the most accurate matching results.
                       For this, we used bulk RNA-seq data from murine and human kidney primary cells as positive controls and matched these to our references.
                       We then assessed how often each method matches our positive control cell type to the correct corresponding reference cell type (figure 1).</p>
                       <p>Spearman's correlation, which compares the rank order of genes in between samples, performed slightly better than Euclidean distance
                       which is based on the sum of pairwise gene-count (CPM) differences between samples.
                       Reference Ransick <em>et al.</em> outperformed all other references.
                       Spearman's correlation using tubular marker gene expression of Ransick <em>et al.</em> against all positive controls matched 89% of the positive controls to the correct cell type.
                       Park <em>et al.</em> performed best using global gene expression and Spearman's correlation.
                       Human snRNA-seq reference by Kidney Precision Medicine Project (KPMP) performed best using Euclidean distance and tubular marker gene expression.
                       16 out of 18 positive controls were derived from mice. KPMP might perform better as a reference dataset with human samples.</p>")
                ),
                column(6,
                  tags$img(src = "score_graph.png", height = "100%", style = "max-width: 100%;"),
                  tags$p("Figure 1: \"Positive controls\" (n = 18) from bulk RNA-seq kidney primary cells and primary tissue were tested across 3 references, 2 algorithms and 3 gene sets.
                         Score of 100% means that CellMatchR always matched the correct reference cell type to the positive control cell type.", style = "font-size: 12px;")
                )
              )
            )
          )
        ),
        fluidRow(
          column(10,
            div(style = "background-color: white; padding: 20px;",
              h3("Legal notice (Impressum)", style = "color: #0072B2;"),
              hr(),
              h4("Responsible for the content of this website:"),
              p("Service Project S1"),
              p("Collaborative Research Center 1453"),
              p("Nephrogenetics (NephGen)"),
              p(tags$a(href = "https://www.sfb1453.uni-freiburg.de/", target = "_blank", "https://www.sfb1453.uni-freiburg.de/")),
              br(),
              p("Institute of Epidemiology and Prevention"),
              p("Faculty of Medicine and Medical Center"),
              p("University of Freiburg"),
              p("79106 Freiburg, Germany"),
              p(tags$a(href = "https://www.uniklinik-freiburg.de/genetische-epidemiologie.html", target = "_blank", "https://www.uniklinik-freiburg.de/genetische-epidemiologie.html")),
              br(),
              HTML("<u>We are happy to receive feedback. For this, please contact:</u>"),
              p("Mona Schoberth: mona.schoberth@uniklinik-freiburg.de"),
              p("Dr. Stefan Haug: stefan.haug@uniklinik-freiburg.de")
            )
          )
        ),
        fluidRow(
          column(10,
            div(style = "background-color: white; padding: 20px;",
              h4("Disclaimer"),
              p("The author assumes no responsibility for the topicality, correctness, completeness or quality of information provided.
                Liability claims against the author which relate to material or immaterial nature caused by the use or misuse of any information provided
                through the use of incorrect or incomplete information are excluded unless the author is not intentional or grossly negligent fault.
                The author reserves the right to change parts of the site without prior notice, add to, delete or cease publication temporarily or permanently.")
            )
          )
        ),
        fluidRow(
          column(10,
            div(style = "background-color: white; padding: 20px;",
              h4("Shinylive Version"),
              p("This version runs entirely in your browser using WebAssembly. No data is uploaded to any server.")
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# Server
# ============================================================================

server <- function(input, output, session) {
  options(shiny.maxRequestSize = 100 * 1024^2)

  rv_sample <- reactiveValues(data = NULL)
  rv_custom_ref <- reactiveValues(data = NULL)
  rv <- reactiveValues(data = NULL, gene_lists = list())

  # Reference data
  reference <- reactive({
    if (!is.null(rv_custom_ref$data)) {
      return(rv_custom_ref$data)
    }

    switch(input$ref,
           "Ransick et al. (mouse, recommended)" = read.csv("data/Ransick_CPM_2025-08-14.csv"),
           "Park et al. (mouse)" = read.csv("data/Park_CPM_2025-08-14.csv"),
           "Lake et al. (human)" = read.csv("data/KPMP_CPM_2025-08-14.csv"),
           "Zhang et al. (human)" = read.csv("data/Zhang_CPM_2025-08-14.csv"))
  })

  # Custom reference upload
  observeEvent(input$custom_ref, {
    req(input$custom_ref)
    rv_custom_ref$data <- read.csv(input$custom_ref$datapath)
  })

  # Sample upload
  observeEvent(input$test, {
    req(input$test)
    rv_sample$data <- read.csv(input$test$datapath)
    updateSelectInput(session, "demo", selected = "--")
  })

  # Reset sample
  observeEvent(input$reset_test, {
    rv_sample$data <- NULL
    shinyjs::reset("test")
  })

  # Sample data
  sample <- reactive({
    req(input$demo)

    if (!is.null(rv_sample$data) && input$demo == "--") {
      return(rv_sample$data)
    }

    switch(input$demo,
           "kidney primary cells - select cell types below" = read.csv("data/CPM_Chen_primary_cells.csv"),
           "HK-2 proximal tubule cell line" = read.csv("data/TPM_HK2_khundmiri.csv"),
           "mIMCD-3 cell line" = read.csv("data/mIMCD_WT_10d_Westermann.csv"),
           NULL)
  })

  # Gene set
  geneset <- reactive({
    if (!is.null(rv$data)) {
      return(rv$data)
    }

    switch(input$genes,
           "all genes" = NULL,
           "kidney marker genes" = read.csv("data/marker_genes_all_2025-07-16.csv"),
           "tubular marker genes" = read.csv("data/marker_genes_tubular_2025-07-16.csv"))
  })

  # Load gene-celltype assignment for heatmap grouping
  gene_celltype <- read.csv("data/mg_all_celltype_assignment.csv")

  # Custom gene list upload
  observeEvent(input$genelist, {
    req(input$genelist)
    rv$data <- read.csv(input$genelist$datapath)
    rv$gene_lists[[input$genelist$name]] <- rv$data
    updateSelectInput(session, "genes",
                      choices = c("all genes", "kidney marker genes", "tubular marker genes", names(rv$gene_lists)))
  })

  # Update cell type pickers
  names_ref <- reactive({
    colnames(reference())[colnames(reference()) != "Gene.names"]
  })

  names_test <- reactive({
    colnames(sample())[colnames(sample()) != "Gene.names"]
  })

  observeEvent(reference(), {
    updatePickerInput(session, "reference_celltype_select",
                      choices = names_ref(), selected = names_ref())
  })

  observeEvent(sample(), {
    updatePickerInput(session, "test_celltype_select",
                      choices = names_test(), selected = names_test())
  })

  # Selected subsets
  reference_select <- reactive({
    req(input$reference_celltype_select)
    valid_cols <- intersect(input$reference_celltype_select, colnames(reference()))
    req(length(valid_cols) > 0)
    reference()[, c("Gene.names", valid_cols)]
  })

  test_select <- reactive({
    if (is.null(sample())) return(NULL)
    req(input$test_celltype_select)
    valid_cols <- intersect(input$test_celltype_select, colnames(sample()))
    req(length(valid_cols) > 0)
    sample()[, c("Gene.names", valid_cols)]
  })

  selected_option <- reactive({ input$counts })

  # Reduce data
  reduced_data <- reactive({
    if (is.null(rv_sample$data) && input$demo == "--") {
      return(NULL)
    } else if (is.null(rv_sample$data) && input$demo != "--") {
      reduce_data(reference_select(), test_select(), geneset())
    } else if (!is.null(rv_sample$data) && input$demo == "--" && selected_option() == "counts") {
      reduce_data_counts(reference_select(), test_select(), geneset())
    } else if (!is.null(rv_sample$data) && input$demo == "--" && selected_option() == "CPM") {
      reduce_data(reference_select(), test_select(), geneset())
    }
  })

  # Results
  results_cor <- reactive({
    if (is.null(rv_sample$data) && input$demo == "--") return(NULL)
    cor_data(reduced_data())
  })

  # results_ED <- reactive({
  #   if (is.null(rv_sample$data) && input$demo == "--") return(NULL)
  #   ED_data(reduced_data())
  # })

  # Legend
  output$legend <- renderTable({
    if (!is.null(rv_custom_ref$data)) {
      return(data.frame(Note = "Custom reference uploaded - no legend available"))
    }
    if (input$ref %in% names(legends)) {
      legends[[input$ref]]
    }
  })

  # Plots
  output$cor_plot <- renderPlot({
    if (input$button > 0) {
      graph(results_cor())
    }
  })

  # output$ED_plot <- renderPlot({
  #   if (input$button > 0) {
  #     graph_ED(results_ED())
  #   }
  # })

  # Heatmap
  heatmap_result <- reactive({
    req(reduced_data(), results_cor())
    genes_to_plot <- input$expression_genes_input
    scale_opt <- input$scale == "scale by row"
    heatmap_fct(reduced_data(), results_cor(), curated_markers_df,
                genesToPlot = genes_to_plot, scale = scale_opt,
                gene_celltype = gene_celltype)
  })

  heatmap_height <- reactive({
    req(heatmap_result())
    calculate_heatmap_height(heatmap_result()$n_genes)
  })

  output$heatmap_ui <- renderUI({
    h <- paste0(heatmap_height(), "px")
    plotOutput("heatmap", width = "100%", height = h)
  })

  output$heatmap <- renderPlot({
    req(heatmap_result())
    print(heatmap_result()$plot)
  })

  observeEvent(input$reset_genes, {
    updateSelectizeInput(session, "expression_genes_input", selected = character(0))
  })

  # Results tabBox
  output$tabBox <- renderUI({
    if (input$button == 0 || (is.null(rv_sample$data) && input$demo == "--")) return(NULL)

    fluidRow(
      tags$style(type = "text/css",
                 ".shiny-output-error { visibility: hidden; }",
                 ".shiny-output-error:before { visibility: hidden; }"),
      div(
        style = "padding: 18px;",
        tabBox(
          title = "Results",
          width = 12,
          tabPanel("Spearman correlation",
            br(),
            p("Results show Spearman's rho ranked from highest to lowest. Higher rho = more similar.
              Dashed lines at rho = 0.8 (black) and 0.6 (grey). Error bars show min/max across replicates."),
            fluidRow(
              column(8,
                plotOutput("cor_plot", width = "100%", height = "800px"),
                downloadButton("DL_cor_csv", "Download results (CSV)")
              )
            )
          ),
          tabPanel("Heatmap",
            fluidRow(
              column(4,
                selectizeInput("expression_genes_input",
                               label = tags$div("Gene search",
                                                actionLink("reset_genes", "Reset", style = "padding-left: 20px;")),
                               choices = sort(rownames(reduced_data()$ref_filt)),
                               options = list(
                                 maxOptions = 100,
                                 placeholder = "Type gene name(s) to search..."
                               ),
                               multiple = TRUE)
              ),
              column(4,
                selectInput("scale", label = "Heatmap options",
                            choices = c("no scaling", "scale by row"))
              )
            ),
            br(),
            p("By default, heatmap displays", strong("40 selected tubular marker genes."),
              "However, all genes of the datatable can be selected."),
            p("Genes are labelled with the respective cell types and appended with a dot and an incremental number in case they are markers for multiple tubular cell types."),
            hr(),
            fluidRow(
              column(12, uiOutput("heatmap_ui"))
            ),
            fluidRow(
              column(4, downloadButton("DL_heatmap", "Download heatmap data (CSV)"))
            )
          )
        )
      )
    )
  })

  # Downloads
  output$DL_cor_csv <- downloadHandler(
    filename = function() paste0("cellmatchr_correlation_", Sys.Date(), ".csv"),
    content = function(file) write.csv(results_cor(), file, row.names = FALSE)
  )

  # output$DL_ED_csv <- downloadHandler(
  #   filename = function() paste0("cellmatchr_euclidean_", Sys.Date(), ".csv"),
  #   content = function(file) write.csv(results_ED(), file, row.names = FALSE)
  # )

  output$DL_heatmap <- downloadHandler(
    filename = function() paste0("cellmatchr_heatmap_", Sys.Date(), ".csv"),
    content = function(file) {
      req(heatmap_result())
      write.csv(heatmap_result()$download_matrix, file, row.names = FALSE)
    }
  )
}

# ============================================================================
# Run App
# ============================================================================

shinyApp(ui = ui, server = server)
