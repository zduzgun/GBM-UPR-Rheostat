# ==============================================================================
# 07_02_neftel_clustering_and_visualization.R
#
# Purpose: To load the initial Neftel Seurat object, perform clustering,
#          run UMAP for visualization, and save the updated object.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2) # FIX: Load ggplot2 for the ggtitle() function
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
input_rds_file <- "07_01_neftel_initial_seurat_object.rds"
output_rds_file <- "07_02_neftel_clustered_seurat_object.rds"
output_plot_file <- "07_02_neftel_umap_clusters.pdf"

# --- 3. Load Seurat Object ---
cat("Step 3: Loading initial Seurat object...\n")
neftel_seurat <- readRDS(input_rds_file)
cat("Seurat object loaded successfully.\n")

# --- 4. Clustering ---
cat("Step 4: Performing cell clustering...\n")
cat("   - Finding neighbors based on PCA...\n")
neftel_seurat <- FindNeighbors(neftel_seurat, dims = 1:30)

cat("   - Finding clusters...\n")
neftel_seurat <- FindClusters(neftel_seurat, resolution = 0.8)

cat("Clustering complete. Found", length(levels(neftel_seurat)), "clusters.\n")

# --- 5. UMAP Visualization ---
cat("Step 5: Running UMAP for visualization...\n")
neftel_seurat <- RunUMAP(neftel_seurat, dims = 1:30)
cat("UMAP calculation complete.\n")

# --- 6. Generate and Save Plot ---
cat("Step 6: Generating and saving UMAP plot...\n")
p <- DimPlot(neftel_seurat, reduction = "umap", label = TRUE, repel = TRUE) +
     ggtitle("UMAP of Neftel et al. GBM Cells")

pdf(output_plot_file, width = 10, height = 8)
print(p)
dev.off()
cat("UMAP plot saved to:", output_plot_file, "\n")

# --- 7. Save the Result ---
cat("Step 7: Saving the clustered Seurat object...\n")
saveRDS(neftel_seurat, file = output_rds_file)
cat("Analysis complete. Final object saved to '", output_rds_file, "'.\n")