#------------------------------------------#
# graph function, Spearman correlation, Shiny
#-----------------------------------------#

graph <- function(results) {
 
  results$min[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA
  results$max[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA
  
  gg <- ggplot(results, aes(x=rho, y= reorder(celltypes, rho), fill = datatype))+
    geom_col(width = 0.65) +
    geom_text(aes(label=sprintf("%.2f", rho)), vjust = 0.5, hjust = -0.2, size = 4.2, color="black") +
    geom_errorbar(aes(xmin = min, xmax = max), width = 0.2, linewidth = 0.5, color = "grey40") +
    geom_vline(xintercept = 0.8, linetype = "dashed", colour = "black") + 
    geom_vline(xintercept = 0.6, linetype = "dashed", colour = "grey60") +
    scale_x_continuous(breaks = seq(0, 1, by = 0.2),
                       expand = expansion(mult = c(0, 0.05))) +
    expand_limits(x = c(0, 1.05 * max(results$rho))) +
    facet_grid(datatype ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Spearman's rho", y = NULL) +
    scale_fill_manual(values = c("#B85C00", "#3B75A3")) +
    theme_minimal(base_size  = 15) +
    theme(
      axis.text.y = element_text(size = 12, color = "black"),
      axis.text.x = element_text(size = 12),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 12),
      panel.grid.major.y = element_blank(), # removes busy grid lines
      panel.grid.minor = element_blank(),
      strip.text.y = element_blank()        # keeps facets clean
    )

  return(gg)
}

