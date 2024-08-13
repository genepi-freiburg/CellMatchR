# cor data, shiny

cor_data <- function(reduce_data) {
  
  ref_test <- cbind(data.frame(reduce_data$ref_filt[, 1:dim(reduce_data$ref_filt)[2]]),
                    data.frame(reduce_data$test_filt[, 1:dim(reduce_data$test_filt)[2]]),
                    as.data.frame(reduce_data$test_median))
  colnames(ref_test)[ncol(ref_test)] <- "test_median"
  colnames(ref_test)[(ncol(reduce_data$ref_filt) + 1) : (ncol(reduce_data$ref_filt) + ncol(reduce_data$test_filt))] <- colnames(reduce_data$test_filt)
  
  rho <- as.data.frame(cor(ref_test, method = "spearman"))
  rho <- rho[-nrow(rho), ]
  
  min <- as.data.frame(apply(rho[, ((dim(reduce_data$ref_filt)[2] + 1) : ((dim(reduce_data$ref_filt)[2] + 1) + dim(reduce_data$test_filt)[2]))], 1, min)) 
  colnames(min) <- c("min")
  
  max <- as.data.frame(apply(rho[, ((dim(reduce_data$ref_filt)[2] + 1) : ((dim(reduce_data$ref_filt)[2] + 1) + dim(reduce_data$test_filt)[2]))], 1, max))
  colnames(max) <- "max"
  
  cor.df <- data.frame(rho = rho$test_median, celltypes = rownames(rho))
  
  results_df <- mutate(cor.df, datatype = ifelse(row_number() <= length(colnames(reduce_data$ref_filt)), "Median of sample(s) vs. reference", "Median of sample(s) vs. sample(s)"))
  
  duplicates <- duplicated(results_df$celltypes)                                   #wenn duplizierte Namen in celltypes, mach Leerzeichen vor 2.Duplikat
  results_df$celltypes[duplicates] <- paste0(" ", results_df$celltypes[duplicates])
  
  results <- cbind(results_df, min, max)
  
  return(results)
}