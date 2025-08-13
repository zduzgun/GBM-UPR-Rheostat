# ==============================================================================
# 07_05_neftel_cluster_characterization.R
#
# Purpose: To identify marker genes for each cluster in the Neftel dataset
#          to understand their biological identity and compare them to the
#          primary analysis.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
# We use the clustered object before any subsampling was done for trajectory.
input_rds_file <- "07_02_neftel_clustered_seurat_object.rds"
output_csv_file <- "07_05_neftel_all_cluster_markers.csv"
output_plot_file <- "07_05_neftel_top10_markers_heatmap.pdf"

# --- 3. Load Clustered Seurat Object ---
cat("Step 3: Loading clustered Seurat object...\n")
neftel_seurat <- readRDS(input_rds_file)
# Set the identity of the cells to the clusters we found
Idents(neftel_seurat) <- "seurat_clusters"
cat("Seurat object loaded. Active identity set to 'seurat_clusters'.\n")
cat("Total cells to analyze:", ncol(neftel_seurat), "\n")


# --- 4. Find Marker Genes for Each Cluster ---
cat("Step 4: Finding marker genes for each cluster. This may take a while...\n")
# FindAllMarkers will compare each cluster to all other clusters.
all_markers <- FindAllMarkers(
  neftel_seurat,
  only.pos = TRUE,          # Only identify positive markers
  min.pct = 0.25,           # Detect gene in at least 25% of cells in the cluster
  logfc.threshold = 0.25    # Minimum log2 fold change
)
cat("FindAllMarkers complete. Found", nrow(all_markers), "total marker genes.\n")

# --- 5. Save the Full Marker List ---
cat("Step 5: Saving the complete list of marker genes to CSV...\n")
write.csv(all_markers, file = output_csv_file, row.names = FALSE)
cat("Marker gene list saved to:", output_csv_file, "\n")


# --- 6. Generate a Heatmap of Top Markers for Visualization ---
cat("Step 6: Generating a heatmap of top 10 markers per cluster...\n")
# Get the top 10 markers for each cluster based on average log2 fold change
top10_markers <- all_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Create the heatmap
p <- DoHeatmap(neftel_seurat, features = top10_markers$gene) +
     NoLegend() +
     ggtitle("Top 10 Marker Genes per Cluster - Neftel Dataset")

pdf(output_plot_file, width = 16, height = 20)
print(p)
dev.off()
cat("Heatmap plot saved to:", output_plot_file, "\n")
cat("Analysis complete.\n")