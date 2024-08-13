#------------------------------------------------------------#
# Download
#-----------------------------------------------------------#

library(writexl)

# Download graph and table from Spearman's correlation

output$DL_cor_plot <- downloadHandler(
  filename = function() {
    paste0("spearmans_correlation_plot.pdf")
  },
  content = function(file) {
    pdf(file, width = 12, height = 12)
    print(graph(results_cor()))
    dev.off()
  })

output$DL_cor_table <- downloadHandler(
  filename = function() {
    paste0("spearmans_correlation_table.xlsx")
  },
  content = function(file) {
    write_xlsx(results_cor(), file)
  }
)

# Download graph and table from Euclidean distance

output$DL_ED_plot <- downloadHandler(
  filename = function() {
    paste0("euclidean_distance_plot.pdf")
  },
  content = function(file) {
    pdf(file, width = 12, height = 12)
    print(graph_ED(results_ED()))
    dev.off()
  })

output$DL_ED_table <- downloadHandler(
  filename = function() {
    paste0("euclidean_distance_table.xlsx")
  },
  content = function(file) {
    write_xlsx(results_ED(), file)
  }
)

# Download Excel table of heatmap dataframe

output$DL3 <- downloadHandler(
  filename = function() {
    paste0("heatmap_datatable.xlsx")
  },
  content = function(file) {
    write_xlsx(results_map()$download_matrix, file)
  }
)
