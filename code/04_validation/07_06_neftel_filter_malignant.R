# ==============================================================================
# 07_06_neftel_filter_malignant.R
#
# Purpose: To filter out non-malignant cell clusters (immune cells,
#          contaminants) from the Neftel dataset to create a "clean"
#          tumor-only Seurat object for downstream analysis.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
# We use the full clustered object as input
input_rds_file <- "07_02_neftel_clustered_seurat_object.rds"
output_rds_file <- "07_06_neftel_malignant_only.rds"
output_plot_file <- "07_06_neftel_filtering_comparison.pdf"

# --- 3. Load Clustered Seurat Object ---
cat("Step 3: Loading full clustered Seurat object...\n")
neftel_seurat <- readRDS(input_rds_file)
cat("Seurat object loaded successfully.\n")
cat("Total cells before filtering:", ncol(neftel_seurat), "\n")

# --- 4. Define Non-Malignant Clusters to Remove ---
# Based on our marker gene analysis (07_05)
cat("Step 4: Defining non-malignant clusters to be removed...\n")
immune_clusters <- c(2, 3, 4, 12, 16, 21) # Microglia/Macrophages and T-cells
contaminant_clusters <- c(28) # Erythroid cells

clusters_to_remove <- c(immune_clusters, contaminant_clusters)
cat("Clusters to remove:", paste(clusters_to_remove, collapse = ", "), "\n")

# --- 5. Generate "Before" Plot ---
cat("Step 5: Generating UMAP plot before filtering...\n")
# Highlight the clusters that will be removed
neftel_seurat$Status <- ifelse(neftel_seurat$seurat_clusters %in% clusters_to_remove, "ToRemove", "Keep")
p_before <- DimPlot(neftel_seurat, group.by = "Status", cols = c("Keep" = "grey", "ToRemove" = "red")) +
            ggtitle("Before Filtering: Non-Malignant Clusters Highlighted")

# --- 6. Filter the Seurat Object ---
cat("Step 6: Subsetting the data to keep only malignant cells...\n")
malignant_seurat <- subset(neftel_seurat, idents = clusters_to_remove, invert = TRUE)
cat("Total cells after filtering:", ncol(malignant_seurat), "\n")
cat("Removed", ncol(neftel_seurat) - ncol(malignant_seurat), "non-malignant cells.\n")

# --- 7. Generate "After" Plot ---
cat("Step 7: Generating UMAP plot after filtering...\n")
p_after <- DimPlot(malignant_seurat, reduction = "umap", label = TRUE) +
           ggtitle("After Filtering: Malignant & Stromal Cells Only")

# --- 8. Save Comparison Plot ---
cat("Step 8: Saving comparison plot...\n")
comparison_plot <- p_before + p_after

pdf(output_plot_file, width = 16, height = 8)
print(comparison_plot)
dev.off()
cat("Comparison plot saved to:", output_plot_file, "\n")

# --- 9. Save the Cleaned Seurat Object ---
cat("Step 9: Saving the cleaned, malignant-only Seurat object...\n")
# It's good practice to re-run the standard analysis steps on the subset
# to ensure the structure is clean and ready for deep analysis.
# However, for now we will save the direct subset.
saveRDS(malignant_seurat, file = output_rds_file)
cat("Cleaned Seurat object saved to:", output_rds_file, "\n")

cat("Analysis complete.\n")