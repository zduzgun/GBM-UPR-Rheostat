#################################################################
# SCRIPT 06_03: UPR SCORE CORRELATION ANALYSIS
# Purpose: To quantitatively test the "UPR Rheostat" hypothesis
#          by calculating the correlation between the three UPR arms.
# Input:   glioblastoma_with_final_scores.rds
# Output:  upr_correlation_results.txt
#          upr_correlation_heatmap.pdf
#################################################################

# Load required libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(corrplot)

# --- Load Data ---
# Load the Seurat object containing all calculated scores
print("Reading Seurat object: glioblastoma_with_final_scores.rds")
seurat_obj <- readRDS("glioblastoma_with_final_scores.rds")
print("Data loading complete.")

# --- Extract UPR Scores ---
# Extract the metadata dataframe from the Seurat object
metadata <- seurat_obj@meta.data

# Select only the UPR score columns for analysis
# Note: Seurat adds a '1' to the end of the score names
upr_scores_df <- metadata %>%
  select(IRE1_Score1, PERK_Score1, ATF6_Score1)

print("UPR scores extracted for analysis.")

# --- Calculate Correlation Matrices ---
# Calculate both Pearson and Spearman correlations
# Pearson looks for linear relationships
# Spearman looks for monotonic relationships (robust to outliers)
print("Calculating Pearson and Spearman correlation matrices...")
pearson_cor_matrix <- cor(upr_scores_df, method = "pearson")
spearman_cor_matrix <- cor(upr_scores_df, method = "spearman")

# --- Perform Statistical Significance Test for Pearson Correlation ---
# corr.test from the 'psych' package is good for this, but to avoid
# extra dependencies, we can do it manually for each pair.
# A simpler approach is using cor.mtest from the 'corrplot' package.
cor_test_results <- cor.mtest(upr_scores_df, conf.level = 0.95)
p_value_matrix <- cor_test_results$p

print("Correlation analysis complete.")

# --- Save Text Output ---
# Save the matrices and p-values to a text file for review
sink("upr_correlation_results.txt")
print("--- Pearson Correlation Matrix ---")
print(pearson_cor_matrix)
print("\n--- Spearman Correlation Matrix ---")
print(spearman_cor_matrix)
print("\n--- P-values for Pearson Correlation ---")
print(p_value_matrix)
sink()

print("Correlation results saved to upr_correlation_results.txt")


# --- Generate Visualization ---
# Create a correlation heatmap to visualize the results for the paper (Figure 2C)
print("Generating correlation heatmap...")
pdf("upr_correlation_heatmap.pdf", width = 8, height = 8)
corrplot(pearson_cor_matrix,
         method = "color",       # Use color to represent correlation
         type = "upper",         # Show only the upper triangle
         order = "hclust",       # Reorder based on hierarchical clustering
         addCoef.col = "black",  # Add correlation coefficients
         tl.col = "black",       # Text label color
         tl.srt = 45,            # Rotate text labels
         p.mat = p_value_matrix, # Provide the p-value matrix
         sig.level = 0.05,       # Significance level
         insig = "blank",        # Leave non-significant correlations blank
         diag = FALSE,           # Don't show the diagonal
         title = "Correlation Matrix of UPR Pathway Scores",
         mar = c(0,0,1,0))       # Adjust margins
dev.off()

print("Heatmap saved to upr_correlation_heatmap.pdf")
print("Process completed successfully.")