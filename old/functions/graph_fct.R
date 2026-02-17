#------------------------------------------#
# graph function, Spearman correlation, Shiny
#-----------------------------------------#

graph <- function(results) {
 
  results$min[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA
  results$max[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA
  
  gg <- ggplot(results, aes(x=rho, y= reorder(celltypes, rho), fill = datatype, width = 0.7))+
    geom_bar(stat="identity") +
    geom_text(aes(label=round(rho, digits = 2)), vjust = 0.5, size = 5.0, hjust = -0.8, color="black") +
    geom_errorbar(aes(xmin = min, xmax = max), width = 0.3, linewidth = 0.75) +
    labs(title = NULL , x = "Spearman's rho") +
    expand_limits(x = c(0, 1.30 * max(results$rho))) +
    theme_light() +
    facet_grid(datatype ~ ., scales = "free_y", space = "free_y") +
    theme(plot.title = element_text(hjust = 0, size = 14),
          axis.title.y = element_blank(),
          strip.text.y = element_blank(),
          text = element_text(size = 15),
          legend.position = "bottom",
          legend.direction = "vertical",
          legend.title=element_blank(),
          legend.text = element_text(size = 11))+
    scale_fill_manual(values = c("#bd4861", "#0072B2")) 
  
  return(gg)
}

