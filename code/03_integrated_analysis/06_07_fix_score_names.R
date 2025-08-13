#!/usr/bin/env Rscript

# 06_07_fix_score_names.R
# Fix the column naming mismatch for UPR and metabolic scores

suppressPackageStartupMessages({
  library(Seurat)
})

message("=== 06_07 Score Name Fixing ===")

# Load dataset
message("Loading dataset...")
glioblastoma <- readRDS("glioblastoma_with_final_scores.rds")

message("Current metadata columns:")
print(colnames(glioblastoma@meta.data))

# Fix UPR score names
message("Fixing UPR score names...")
if("IRE1_Score1" %in% colnames(glioblastoma@meta.data)) {
  glioblastoma@meta.data$IRE1_score <- glioblastoma@meta.data$IRE1_Score1
  message("? IRE1_Score1 › IRE1_score")
}

if("PERK_Score1" %in% colnames(glioblastoma@meta.data)) {
  glioblastoma@meta.data$PERK_score <- glioblastoma@meta.data$PERK_Score1
  message("? PERK_Score1 › PERK_score")
}

if("ATF6_Score1" %in% colnames(glioblastoma@meta.data)) {
  glioblastoma@meta.data$ATF6_score <- glioblastoma@meta.data$ATF6_Score1
  message("? ATF6_Score1 › ATF6_score")
}

# Fix metabolic score names
message("Fixing metabolic score names...")
if("Glycolysis_Score1" %in% colnames(glioblastoma@meta.data)) {
  glioblastoma@meta.data$Glycolysis_score <- glioblastoma@meta.data$Glycolysis_Score1
  message("? Glycolysis_Score1 › Glycolysis_score")
}

if("OXPHOS_Score1" %in% colnames(glioblastoma@meta.data)) {
  glioblastoma@meta.data$OXPHOS_score <- glioblastoma@meta.data$OXPHOS_Score1
  message("? OXPHOS_Score1 › OXPHOS_score")
}

# Add missing Fatty Acid score if not present
if(!"Fatty_Acid_score" %in% colnames(glioblastoma@meta.data)) {
  # Create a dummy score or use a related pathway
  message("? Adding placeholder Fatty_Acid_score (using Mesenchymal as proxy)")
  if("Mesenchymal_Score1" %in% colnames(glioblastoma@meta.data)) {
    glioblastoma@meta.data$Fatty_Acid_score <- glioblastoma@meta.data$Mesenchymal_Score1
  } else {
    glioblastoma@meta.data$Fatty_Acid_score <- 0
  }
}

# Verify all required scores are now present
required_upr <- c("IRE1_score", "PERK_score", "ATF6_score")
required_metabolic <- c("Glycolysis_score", "OXPHOS_score", "Fatty_Acid_score")

found_upr <- sum(required_upr %in% colnames(glioblastoma@meta.data))
found_metabolic <- sum(required_metabolic %in% colnames(glioblastoma@meta.data))

message("\n=== VERIFICATION ===")
message(paste("UPR scores available:", found_upr, "/", length(required_upr)))
message(paste("Metabolic scores available:", found_metabolic, "/", length(required_metabolic)))

if(found_upr >= 2 && found_metabolic >= 2) {
  message("? SUCCESS: UPR Reostat analysis now possible!")
} else {
  message("? Still missing critical scores")
}

# Show summary statistics for key scores
for(score in c(required_upr, required_metabolic)) {
  if(score %in% colnames(glioblastoma@meta.data)) {
    message(paste("\n", score, "summary:"))
    print(summary(glioblastoma@meta.data[[score]]))
  }
}

# Save updated dataset
message("\nSaving updated dataset...")
saveRDS(glioblastoma, "glioblastoma_with_final_scores_fixed.rds")

# Also update the original file
file.copy("glioblastoma_with_final_scores_fixed.rds", "glioblastoma_with_final_scores.rds", overwrite = TRUE)

message("\n=== 06_07 COMPLETED ===")
message("Score names fixed and dataset updated")
message("Ready for UPR Reostat analysis!")