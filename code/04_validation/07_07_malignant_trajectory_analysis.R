# ==============================================================================
# 07_07_malignant_trajectory_analysis.R (v2 - Monocle Clustering Fix)
#
# Purpose: Perform a refined trajectory analysis using only the malignant
#          cell subset from the Neftel dataset. This version fixes the
#          missing cluster_cells step required by learn_graph.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
input_rds_file <- "07_06_neftel_malignant_only.rds"
output_cds_file <- "07_07_malignant_monocle_cds.rds"
output_plot_file <- "07_07_malignant_trajectory_plots.pdf"
output_seurat_rds_file <- "07_07_malignant_processed_seurat.rds"

# --- 3. Load Filtered Seurat Object ---
cat("Step 3: Loading malignant-only Seurat object...\n")
seurat_obj <- readRDS(input_rds_file)
cat("Seurat object loaded. Contains", ncol(seurat_obj), "malignant cells.\n")

# --- 4. Re-run Seurat Workflow on Subset ---
cat("Step 4: Re-running PCA and UMAP on the malignant-only subset...\n")
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)
seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
cat("Re-processing complete.\n")

# --- 5. Convert Seurat to Monocle3 CDS object ---
cat("Step 5: Converting re-processed Seurat object to Monocle3 format...\n")
expression_matrix <- GetAssayData(seurat_obj, assay = 'RNA', layer = 'counts')
cell_metadata <- seurat_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), row.names = rownames(expression_matrix))

cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

# --- 6. Process CDS and Learn Trajectory ---
cat("Step 6: Preprocessing CDS and learning trajectory...\n")
# Transfer the NEW PCA and UMAP
reducedDims(cds)[['PCA']] <- seurat_obj@reductions$pca@cell.embeddings
reducedDims(cds)[['UMAP']] <- seurat_obj@reductions$umap@cell.embeddings

# ================================ FIX IS HERE =================================
# CRITICAL FIX: Run cluster_cells first, then optionally override with Seurat clusters
cat("  - Running Monocle3 clustering...\n")
cds <- cluster_cells(cds, reduction_method = "UMAP")

# Optional: If you want to use Seurat clusters instead of Monocle3 clusters
cat("  - Optionally overriding with Seurat clusters...\n")
colData(cds)$seurat_clusters <- seurat_obj@meta.data$seurat_clusters
colData(cds)$assigned_cluster <- seurat_obj@meta.data$seurat_clusters
# ===============================================================================

# FIXED: Learn the graph with more flexible parameters  
cat("  - Learning the trajectory graph with flexible parameters...\n")
cds <- learn_graph(cds, 
                   use_partition = FALSE,           # Don't respect partitions - force connections
                   close_loop = FALSE,              # Don't close loops
                   learn_graph_control = list(
                     minimal_branch_len = 5,           # Allow shorter branches (default: 10)
                     euclidean_distance_ratio = 2,     # More flexible distance (default: 1)  
                     geodesic_distance_ratio = 1/3     # Keep default geodesic ratio
                   ))
cat("Graph learned successfully.\n")

# --- 7. Order Cells (Define Pseudotime) ---
cat("Step 7: Ordering cells and defining pseudotime root...\n")
stemness_genes <- c('MYCN', 'MAGEA3', 'MAGEA12', 'SOX2', 'PROM1')
genes_in_data <- intersect(stemness_genes, rownames(cds))

if (length(genes_in_data) > 0) {
  # AddModuleScore must be run on the Seurat object
  seurat_obj <- AddModuleScore(seurat_obj, features = list(genes_in_data), name = "Stemness_Score")
  colData(cds)$Stemness_Score <- seurat_obj$Stemness_Score1

  # Use the actual Monocle3 clusters, not the Seurat clusters
  colData(cds)$monocle_clusters <- clusters(cds)
  
  root_cluster_info <- as.data.frame(colData(cds)) %>%
    group_by(monocle_clusters) %>%
    summarise(mean_stemness = mean(Stemness_Score, na.rm = TRUE)) %>%
    filter(!is.na(monocle_clusters)) %>%
    arrange(desc(mean_stemness))
  
  root_cluster <- root_cluster_info$monocle_clusters[1]
  cat("  - Identified root cluster based on max stemness score:", root_cluster, "\n")
  
  # Find root cells using the correct Monocle3 cluster IDs - use cell names, not indices
  root_cell_indices <- which(colData(cds)[, "monocle_clusters"] == root_cluster)
  root_cells <- rownames(colData(cds))[root_cell_indices]
  cat("  - Found", length(root_cells), "root cells in cluster", root_cluster, "\n")

  cds <- order_cells(cds, root_cells = root_cells)
  colData(cds)$pseudotime <- pseudotime(cds)
  cat("Pseudotime calculated successfully.\n")
} else {
  cat("  - WARNING: No stemness marker genes found. Letting Monocle choose root.\n")
  cds <- order_cells(cds)
  colData(cds)$pseudotime <- pseudotime(cds)
}

# --- 8. Generate and Save Plots ---
cat("Step 8: Generating and saving refined trajectory plots...\n")
pdf(output_plot_file, width = 12, height = 10)

p_clusters <- plot_cells(cds, color_cells_by = "assigned_cluster", label_groups_by_cluster = TRUE,
                         label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.7) +
              ggtitle("Malignant-Only Trajectory by Cluster")
print(p_clusters)

p_monocle_clusters <- plot_cells(cds, color_cells_by = "monocle_clusters", label_groups_by_cluster = TRUE,
                                 label_leaves = FALSE, label_branch_points = FALSE, cell_size = 0.7) +
                      ggtitle("Malignant-Only Trajectory by Monocle3 Clusters")
print(p_monocle_clusters)

p_pseudotime <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE,
                           label_leaves = TRUE, label_branch_points = TRUE, cell_size = 0.7) +
                ggtitle("Malignant-Only Trajectory by Pseudotime")
print(p_pseudotime)

# Bu satýrlarý kaldýrýn veya yorum satýrý yapýn:
# p_ire1 <- plot_cells(cds, color_cells_by = "IRE1_score1", cell_size = 0.7) +
#           ggtitle("Malignant-Only Trajectory by IRE1 Score") + scale_color_viridis_c()
# print(p_ire1)

# p_perk <- plot_cells(cds, color_cells_by = "PERK_score1", cell_size = 0.7) +
#           ggtitle("Malignant-Only Trajectory by PERK Score") + scale_color_viridis_c()
# print(p_perk)

# p_atf6 <- plot_cells(cds, color_cells_by = "ATF6_score1", cell_size = 0.7) +
#           ggtitle("Malignant-Only Trajectory by ATF6 Score") + scale_color_viridis_c()
# print(p_atf6)

# Bunun yerine Stemness Score plotunu ekleyin:
p_stemness <- plot_cells(cds, color_cells_by = "Stemness_Score", cell_size = 0.7) +
              ggtitle("Malignant-Only Trajectory by Stemness Score") + scale_color_viridis_c()
print(p_stemness)

dev.off()
cat("Plots saved to:", output_plot_file, "\n")

# --- 9. Save Final Objects ---
cat("Step 9: Saving the final Monocle and Seurat objects...\n")
saveRDS(cds, file = output_cds_file)
saveRDS(seurat_obj, file = output_seurat_rds_file)
cat("Analysis complete. Final objects saved successfully.\n")