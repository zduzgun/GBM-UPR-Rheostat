#!/usr/bin/env Rscript

# UPR skorlarý kontrol scripti
suppressPackageStartupMessages({
  library(Seurat)
})

# Ana veri setini yükle
message("Loading glioblastoma dataset...")
glioblastoma <- readRDS("glioblastoma_with_final_scores.rds")

message("Dataset dimensions:", paste(dim(glioblastoma), collapse = " x "))
message("Metadata columns:")
print(colnames(glioblastoma@meta.data))

# UPR skorlarý var mý kontrol et
upr_cols <- c("IRE1_score", "PERK_score", "ATF6_score")
found_upr <- intersect(upr_cols, colnames(glioblastoma@meta.data))

message("Expected UPR scores:", paste(upr_cols, collapse = ", "))
message("Found UPR scores:", paste(found_upr, collapse = ", "))

if(length(found_upr) == 0) {
  message("? NO UPR SCORES FOUND!")
  message("06_01_module_scores.R probably failed or didn't save properly")
} else {
  message("? Found UPR scores:", length(found_upr), "out of", length(upr_cols))
  
  # UPR skorlarý varsa özet istatistik
  for(score in found_upr) {
    message(paste(score, "summary:"))
    print(summary(glioblastoma@meta.data[[score]]))
  }
}

# Metabolik skorlarý da kontrol et
metabolic_cols <- c("Glycolysis_score", "OXPHOS_score", "Fatty_Acid_score")
found_metabolic <- intersect(metabolic_cols, colnames(glioblastoma@meta.data))

message("Expected metabolic scores:", paste(metabolic_cols, collapse = ", "))
message("Found metabolic scores:", paste(found_metabolic, collapse = ", "))

if(length(found_metabolic) == 0) {
  message("? NO METABOLIC SCORES FOUND!")
} else {
  message("? Found metabolic scores:", length(found_metabolic))
}

# Genel durum özeti
message("\n=== OVERALL STATUS ===")
message("UPR Reostat Analysis Possible:", ifelse(length(found_upr) >= 2, "? YES", "? NO"))
message("Metabolic Integration Possible:", ifelse(length(found_metabolic) >= 1, "? YES", "? NO"))

if(length(found_upr) < 2) {
  message("\n??  CRITICAL: Need to re-run 06_01_module_scores.R")
  message("The UPR reostat hypothesis CANNOT be tested without UPR scores!")
}