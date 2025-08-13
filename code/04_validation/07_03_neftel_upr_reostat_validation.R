# ==============================================================================
# 07_03_neftel_upr_reostat_validation.R (v2 - Syntax Fix)
#
# Purpose: Validate the UPR Reostat hypothesis on the Neftel et al. dataset.
#          This version uses a more robust way of handling column names
#          to prevent syntax errors.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
input_rds_file <- "07_02_neftel_clustered_seurat_object.rds"
output_rds_file <- "07_03_neftel_with_scores.rds"
output_corr_file <- "07_03_neftel_upr_correlation_results.txt"
output_plot_file <- "07_03_neftel_upr_scores_umap.pdf"

# --- 3. Load Clustered Seurat Object ---
cat("Step 3: Loading clustered Seurat object...\n")
neftel_seurat <- readRDS(input_rds_file)
cat("Seurat object loaded successfully.\n")

# --- 4. Define UPR Gene Sets ---
cat("Step 4: Defining UPR gene sets...\n")
upr_gene_sets <- list(
  IRE1_score = c("XBP1", "DNAJB9", "EDEM1", "HSPA5", "PDIA4", "ERO1A", "ERN1"),
  PERK_score = c("ATF4", "DDIT3", "PPP1R15A", "ASNS", "TRIB3", "EIF2AK3"),
  ATF6_score = c("ATF6", "HSPA5", "HSP90B1", "CALR")
)

# --- 5. Calculate Module Scores ---
cat("Step 5: Calculating UPR module scores...\n")
# We will use the default names Seurat creates (e.g., IRE1_score1)
for (set_name in names(upr_gene_sets)) {
  cat("  - Calculating", set_name, "\n")
  genes <- upr_gene_sets[[set_name]]
  genes_in_data <- intersect(genes, rownames(neftel_seurat))
  
  if (length(genes_in_data) == 0) {
      cat("  - WARNING: No genes from", set_name, "found. Skipping.\n")
      # Add a dummy column with the expected name to prevent errors
      neftel_seurat@meta.data[[paste0(set_name, '1')]] <- 0
      next
  }
  
  neftel_seurat <- AddModuleScore(
    object = neftel_seurat,
    features = list(genes_in_data),
    name = set_name,
    ctrl = 100
  )
}
cat("Module scores calculated.\n")
cat("Checking metadata column names after scoring:\n")
print(head(colnames(neftel_seurat@meta.data)))

# --- 6. Visualize Scores on UMAP ---
cat("Step 6: Generating UMAP plots for UPR scores...\n")
# Define the actual column names created by Seurat
score_cols_to_plot <- c("IRE1_score1", "PERK_score1", "ATF6_score1")
upr_plots <- list()

for (score in score_cols_to_plot) {
  p <- FeaturePlot(neftel_seurat, features = score, pt.size = 0.1) +
       ggtitle(paste("Neftel Dataset -", score))
  upr_plots[[score]] <- p
}

pdf(output_plot_file, width = 18, height = 6)
print(plot_grid(plotlist = upr_plots, ncol = 3))
dev.off()
cat("UMAP plots saved to:", output_plot_file, "\n")


# --- 7. Perform Correlation Analysis ---
cat("Step 7: Performing correlation analysis for UPR arms...\n")
# Use the actual column names created by Seurat
score_cols_for_corr <- c("IRE1_score1", "PERK_score1", "ATF6_score1")
upr_scores_df <- neftel_seurat@meta.data[, score_cols_for_corr]

corr_matrix <- cor(upr_scores_df, method = "pearson")

# Critical validation test using the actual column names
corr_test <- cor.test(upr_scores_df$PERK_score1, upr_scores_df$IRE1_score1, method = "pearson")

cat("Correlation analysis complete.\n")

# --- 8. Save Correlation Results ---
cat("Step 8: Saving correlation results to text file...\n")
sink(output_corr_file)
cat("====================================================\n")
cat("UPR Reostat Correlation Validation - Neftel Dataset\n")
cat("====================================================\n\n")
cat("Pearson Correlation Matrix:\n")
print(corr_matrix)
cat("\n\n----------------------------------------------------\n\n")
cat("Detailed Test for PERK vs. IRE1 Correlation:\n")
print(corr_test)
cat("\n\n")
cat("This result will be compared to the primary finding (r approx -0.187)\n")
sink()
cat("Results saved to:", output_corr_file, "\n")

# --- 9. Save Final Seurat Object ---
cat("Step 9: Saving the final Seurat object with scores...\n")
saveRDS(neftel_seurat, file = output_rds_file)
cat("Analysis complete. Final object saved to '", output_rds_file, "'.\n")