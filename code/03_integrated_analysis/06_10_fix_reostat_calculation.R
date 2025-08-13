#!/usr/bin/env Rscript

# 06_10_fix_reostat_calculation.R
# Fix UPR Reostat Index calculation and create corrected figures

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

message("=== Fixing UPR Reostat Index Calculation ===")

# Load data
glioblastoma <- readRDS("glioblastoma_with_final_scores.rds")

# Check UPR score ranges
upr_scores <- c("IRE1_score", "PERK_score", "ATF6_score")
for(score in upr_scores) {
  message(paste(score, "range:", 
                round(min(glioblastoma@meta.data[[score]], na.rm = TRUE), 3), 
                "to", 
                round(max(glioblastoma@meta.data[[score]], na.rm = TRUE), 3)))
}

# CORRECT Reostat Index Calculation
upr_matrix <- as.matrix(glioblastoma@meta.data[, upr_scores])

# Method 1: Standard deviation (most common)
reostat_sd <- apply(upr_matrix, 1, function(x) sd(x, na.rm = TRUE))

# Method 2: Range (max - min)
reostat_range <- apply(upr_matrix, 1, function(x) max(x, na.rm = TRUE) - min(x, na.rm = TRUE))

# Method 3: Coefficient of variation (CV) - CORRECTED
# CV only for positive values, use absolute mean
reostat_cv_corrected <- apply(upr_matrix, 1, function(x) {
  x_shifted <- x - min(x, na.rm = TRUE) + 0.01  # Shift to positive
  sd(x_shifted, na.rm = TRUE) / mean(x_shifted, na.rm = TRUE)
})

# Add to metadata
glioblastoma$UPR_Reostat_SD <- reostat_sd
glioblastoma$UPR_Reostat_Range <- reostat_range  
glioblastoma$UPR_Reostat_CV_Fixed <- reostat_cv_corrected

# Check ranges of corrected indices
message("\nCorrected Reostat Index Ranges:")
message(paste("SD:", round(range(reostat_sd, na.rm = TRUE), 3)))
message(paste("Range:", round(range(reostat_range, na.rm = TRUE), 3)))
message(paste("CV Fixed:", round(range(reostat_cv_corrected, na.rm = TRUE), 3)))

# Summary by cluster
cluster_summary <- glioblastoma@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    IRE1_mean = mean(IRE1_score, na.rm = TRUE),
    PERK_mean = mean(PERK_score, na.rm = TRUE), 
    ATF6_mean = mean(ATF6_score, na.rm = TRUE),
    Reostat_SD = mean(UPR_Reostat_SD, na.rm = TRUE),
    Reostat_Range = mean(UPR_Reostat_Range, na.rm = TRUE),
    Reostat_CV = mean(UPR_Reostat_CV_Fixed, na.rm = TRUE),
    .groups = 'drop'
  )

print("Cluster-wise UPR Summary:")
print(cluster_summary)

# Save corrected data
saveRDS(glioblastoma, "glioblastoma_with_corrected_reostat.rds")
write.csv(cluster_summary, "UPR_cluster_corrected_summary.csv", row.names = FALSE)

# Quick corrected visualization
pdf("Corrected_UPR_Reostat_Index.pdf", width = 12, height = 8)

# Panel A: Corrected Reostat Index (using SD)
p1 <- FeaturePlot(glioblastoma, features = "UPR_Reostat_SD", 
                  pt.size = 0.3, cols = c("lightblue", "red")) +
  labs(title = "Corrected UPR Reostat Index (Standard Deviation)",
       subtitle = "Pathway variability within cells") +
  theme_minimal()

# Panel B: Reostat by cluster
p2 <- ggplot(cluster_summary, aes(x = factor(seurat_clusters), y = Reostat_SD,
                                 fill = factor(seurat_clusters))) +
  geom_col(alpha = 0.8, width = 0.7) +
  scale_fill_brewer(type = "qual", palette = "Set1") +
  labs(title = "UPR Reostat Index by Cluster (Corrected)",
       x = "Cell Cluster", y = "UPR Reostat Index (SD)") +
  theme_minimal() +
  guides(fill = "none")

# Panel C: UPR pathway balance visualization  
upr_long <- glioblastoma@meta.data %>%
  select(seurat_clusters, IRE1_score, PERK_score, ATF6_score) %>%
  reshape2::melt(id.vars = "seurat_clusters", 
                variable.name = "UPR_Pathway", value.name = "Activity")

p3 <- ggplot(upr_long, aes(x = factor(seurat_clusters), y = Activity, 
                          fill = UPR_Pathway)) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = c("IRE1_score" = "#E31A1C", 
                              "PERK_score" = "#1F78B4", 
                              "ATF6_score" = "#33A02C"),
                   labels = c("IRE1?", "PERK", "ATF6")) +
  labs(title = "UPR Pathway Balance Across Clusters",
       x = "Cell Cluster", y = "UPR Activity Score") +
  theme_minimal()

# Combine plots
combined <- p1 / (p2 + p3)
print(combined)

dev.off()

message("\n=== Analysis of Reostat Support ===")

# Test if reostat model is supported
correlations <- cor(glioblastoma@meta.data[, upr_scores], use = "complete.obs")
mean_cor <- mean(correlations[upper.tri(correlations)])

message(paste("Mean UPR correlation:", round(mean_cor, 3)))
message("Reostat interpretation:")
if(mean_cor > 0.7) {
  message("- HIGH correlation (>0.7): UPR pathways co-activate (NOT reostat)")
} else if(mean_cor > 0.3) {
  message("- MODERATE correlation (0.3-0.7): Partial co-regulation (SUPPORTS reostat)")
} else {
  message("- LOW correlation (<0.3): Independent regulation (STRONG reostat)")
}

# Check pathway dominance patterns
pathway_dominance <- glioblastoma@meta.data %>%
  rowwise() %>%
  mutate(
    max_pathway = which.max(c(IRE1_score, PERK_score, ATF6_score)),
    max_pathway_name = c("IRE1", "PERK", "ATF6")[max_pathway]
  ) %>%
  group_by(seurat_clusters, max_pathway_name) %>%
  summarise(n = n(), .groups = 'drop') %>%
  group_by(seurat_clusters) %>%
  mutate(prop = n / sum(n))

message("\nPathway dominance patterns by cluster:")
dominant_patterns <- pathway_dominance %>%
  filter(prop > 0.4) %>%  # >40% dominance
  arrange(seurat_clusters)

if(nrow(dominant_patterns) > 0) {
  print(dominant_patterns)
  message("CONCLUSION: Different clusters show pathway preferences - SUPPORTS reostat model")
} else {
  message("CONCLUSION: No clear pathway dominance patterns - LIMITED reostat support")
}

message("\n=== CORRECTED REOSTAT ANALYSIS COMPLETED ===")
message("Files generated:")
message("- glioblastoma_with_corrected_reostat.rds")
message("- UPR_cluster_corrected_summary.csv") 
message("- Corrected_UPR_Reostat_Index.pdf")