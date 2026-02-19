#!/usr/bin/env Rscript
# Script to convert Excel legends to CSV for Shinylive app
# Run this before each Shinylive deployment to sync legends from original app

library(openxlsx)

# Paths (run from repo root: 13_CellMatchR_App/)
source_excel <- "data/Abbreviations_Celltypes.xlsx"
target_dir <- "CellMatchRLive/data"

# Create target directory if needed
if (!dir.exists(target_dir)) {
  dir.create(target_dir, recursive = TRUE)
}

# Get all sheets
sheets <- getSheetNames(source_excel)

cat("Converting Excel legends to CSV...\n")

for (sheet in sheets) {
  df <- read.xlsx(source_excel, sheet = sheet)

  # Create clean filename
  fname <- gsub(" ", "_", sheet)
  fname <- gsub("\\.", "", fname)
  output_file <- file.path(target_dir, paste0("legend_", fname, ".csv"))

  write.csv(df, output_file, row.names = FALSE)
  cat("  Created:", output_file, "(", nrow(df), "rows )\n")
}

cat("\nDone! Legend CSVs are ready in", target_dir, "\n")
