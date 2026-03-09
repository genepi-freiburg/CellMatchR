#-------------------------------------------------------#
## Euclidean Distance Function - Shiny
#-------------------------------------------------------#

ED_data <- function(reduce_data) {
  ref_test <- rbind(t(log(reduce_data$ref_filt[, 1:dim(reduce_data$ref_filt)[2]] + 1)),
                    t(log(reduce_data$test_filt[, 1:dim(reduce_data$test_filt)[2]] + 1)),
                    t(log(as.data.frame(reduce_data$test_median) + 1)))
  
  dis <- dist(ref_test, method = "euclidean")
  
  dis <- as.data.frame(as.matrix(dis))
  colnames(dis)[ncol(dis)] <- "test_median"
  rownames(dis)[nrow(dis)] <- "test_median"
  dis <- dis[-nrow(dis), ]
  
  med <- as.data.frame(dis[1:(dim(reduce_data$ref_filt)[2] + dim(reduce_data$test_filt)[2]), (dim(reduce_data$ref_filt)[2]) + (dim(reduce_data$test_filt)[2] + 1), drop = FALSE])
  colnames(med)[1] <- "distance"
  med$celltypes <- rownames(med)
  rownames(med) <- NULL
  
  min <- as.data.frame(apply(dis[, (ncol(reduce_data$ref_filt) + 1):(ncol(reduce_data$ref_filt)) + (ncol(reduce_data$test_filt))], 1, min))
  colnames(min) <- c("min")
  
  max <- as.data.frame(apply(dis[, (ncol(reduce_data$ref_filt) + 1):(ncol(reduce_data$ref_filt)) + (ncol(reduce_data$test_filt))], 1, max))
  colnames(max) <- c("max")
  
  results_df <- mutate(med, datatype = ifelse(row_number()  <= length(colnames(reduce_data$ref_filt)), "Median of sample(s) vs. reference", "Median of sample(s) vs. sample(s)"))
  
  duplicates <- duplicated(results_df$celltypes)                                   #wenn duplizierte Namen in celltypes, mach Leerzeichen vor 2.Duplikat
  results_df$celltypes[duplicates] <- paste0(" ", results_df$celltypes[duplicates])
  
  results <- cbind(results_df, min, max)
  
  return(results)
}

graph_ED <- function(results) {
  
  results$min[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA
  results$max[results$datatype == "Median of sample(s) vs. sample(s)"] <- NA
  
  gg <- ggplot(results, aes(x=results$distance, y= reorder(results$celltypes, -results$distance), fill = results$datatype))+
    geom_col(width = 0.65) +
    geom_text(aes(label=sprintf("%.2f", results$distance)), vjust = 0.5, hjust = -0.2, size = 4.2, color="black") +
    geom_errorbar(aes(xmin = min, xmax = max), width = 0.2, linewidth = 0.5, color = "grey40") +
    expand_limits(x = c(0, 1.15 * max(results$distance))) +
    facet_grid(datatype ~ ., scales = "free_y", space = "free_y") +
    labs(x = "Euclidian distance", y = NULL) +
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
