#!/usr/bin/env Rscript

# LIANA Debug Script - Check what's available
suppressPackageStartupMessages({
  library(Seurat)
  library(liana)
  library(dplyr)
})

# Load dataset
message("Loading dataset...")
glioblastoma <- readRDS("glioblastoma_with_final_scores.rds")

# Check LIANA installation and available options
message("\n=== LIANA Debug Information ===")
message(paste("LIANA version:", packageVersion("liana")))

message("\nAvailable resources:")
resources <- show_resources()
print(resources)

message("\nAvailable methods:")  
methods <- show_methods()
print(methods)

# Test simple LIANA call with first available options
message("\n=== Testing Basic LIANA Call ===")

# Use the simplest possible call
test_result <- tryCatch({
  liana_wrap(
    glioblastoma,
    method = methods[1],       # First available method
    resource = resources[1],   # First available resource
    idents_col = 'seurat_clusters',
    verbose = TRUE
  )
}, error = function(e) {
  message(paste("Error:", e$message))
  return(NULL)
})

if (!is.null(test_result)) {
  message("SUCCESS! Basic LIANA call worked.")
  message(paste("Results dimensions:", paste(dim(test_result), collapse = " x ")))
  message("Column names:")
  print(colnames(test_result))
  
  # Save test results
  write.csv(test_result, "liana_debug_results.csv", row.names = FALSE)
  message("Test results saved to: liana_debug_results.csv")
  
} else {
  message("FAILED! Could not run basic LIANA analysis.")
  message("Please check LIANA installation:")
  message("devtools::install_github('saezlab/liana')")
}

message("\n=== Debug completed ===")