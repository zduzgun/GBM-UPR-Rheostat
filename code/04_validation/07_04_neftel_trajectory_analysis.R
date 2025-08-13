# ==============================================================================
# 07_04_neftel_trajectory_analysis_fixed.R
#
# Purpose: Perform trajectory analysis on the Neftel dataset using Monocle3
#          with fixes for large dataset handling and improved error management.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
input_rds_file <- "07_03_neftel_with_scores.rds"
output_cds_file <- "07_04_neftel_monocle_cds.rds"
output_plot_file <- "07_04_neftel_trajectory_plots.pdf"
output_log_file <- "07_04_neftel_trajectory_analysis.log"

# Start logging
sink(file = output_log_file, append = FALSE, split = TRUE)

# --- 3. Load Seurat Object ---
cat("Step 3: Loading Seurat object with scores...\n")
seurat_obj <- readRDS(input_rds_file)
cat("Seurat object loaded successfully.\n")
cat("Original dataset dimensions:", dim(seurat_obj), "\n")
cat("Number of cells:", ncol(seurat_obj), "\n")
cat("Number of genes:", nrow(seurat_obj), "\n")

# --- 4. Subsample Data if Too Large ---
cat("Step 4: Checking dataset size and subsampling if necessary...\n")
max_cells <- 8000  # Safe limit for Monocle3
if (ncol(seurat_obj) > max_cells) {
  cat("Dataset is large (", ncol(seurat_obj), " cells). Subsampling to", max_cells, "cells...\n")
  
  # Stratified sampling to maintain cluster proportions
  set.seed(42)
  cell_clusters <- seurat_obj$seurat_clusters
  cells_per_cluster <- table(cell_clusters)
  
  # Calculate sampling proportion
  sampling_prop <- max_cells / ncol(seurat_obj)
  
  # Sample cells from each cluster proportionally
  sampled_cells <- c()
  for (cluster in names(cells_per_cluster)) {
    cluster_cells <- names(cell_clusters[cell_clusters == cluster])
    n_sample <- max(1, round(length(cluster_cells) * sampling_prop))
    n_sample <- min(n_sample, length(cluster_cells))  # Don't exceed available cells
    sampled_cells <- c(sampled_cells, sample(cluster_cells, n_sample))
  }
  
  # Subset the Seurat object
  seurat_obj <- subset(seurat_obj, cells = sampled_cells)
  cat("Subsampled dataset dimensions:", dim(seurat_obj), "\n")
  cat("Final number of cells:", ncol(seurat_obj), "\n")
} else {
  cat("Dataset size is acceptable (", ncol(seurat_obj), " cells). No subsampling needed.\n")
}

# --- 5. Convert Seurat to Monocle3 CDS object ---
cat("Step 5: Converting Seurat object to Monocle3 format...\n")
# Expression data - use RNA assay with counts layer
expression_matrix <- GetAssayData(seurat_obj, assay = 'RNA', layer = 'counts')
# Cell metadata
cell_metadata <- seurat_obj@meta.data
# Gene metadata
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), 
                           row.names = rownames(expression_matrix))

# Create CDS object
cds <- new_cell_data_set(
  expression_matrix,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cat("CDS object created successfully.\n")

# --- 6. Preprocess Data ---
cat("Step 6: Preprocessing data...\n")
# Transfer existing dimensionality reductions from Seurat
if ("pca" %in% names(seurat_obj@reductions)) {
  reducedDims(cds)[['PCA']] <- seurat_obj@reductions$pca@cell.embeddings
  cat("Transferred PCA from Seurat.\n")
}

if ("umap" %in% names(seurat_obj@reductions)) {
  reducedDims(cds)[['UMAP']] <- seurat_obj@reductions$umap@cell.embeddings
  cat("Transferred UMAP from Seurat.\n")
}

# --- 7. Cluster Cells ---
cat("Step 7: Clustering cells...\n")
# Use existing Seurat clusters for consistency
cds@clusters$UMAP$clusters <- seurat_obj$seurat_clusters
colData(cds)$assigned_cluster <- seurat_obj$seurat_clusters

# Create partitions (important for trajectory learning)
cds@clusters$UMAP$partitions <- rep(1, ncol(cds))
names(cds@clusters$UMAP$partitions) <- colnames(cds)

cat("Clustering completed. Number of clusters:", length(unique(seurat_obj$seurat_clusters)), "\n")

# --- 8. Learn Trajectory Graph ---
cat("Step 8: Learning trajectory graph...\n")
tryCatch({
  cds <- learn_graph(cds, use_partition = FALSE)  # Set to FALSE to avoid partition issues
  cat("Graph learned successfully.\n")
}, error = function(e) {
  cat("Error learning graph:", e$message, "\n")
  cat("Trying with different parameters...\n")
  cds <- learn_graph(cds, use_partition = FALSE, 
                     close_loop = FALSE,
                     learn_graph_control = list(minimal_branch_len = 10))
  cat("Graph learned with alternative parameters.\n")
})

# --- 9. Calculate Stemness Score and Define Root ---
cat("Step 9: Calculating stemness scores and defining trajectory root...\n")
stemness_genes <- c('MYCN', 'MAGEA3', 'MAGEA12', 'SOX2', 'PROM1')
genes_in_data <- intersect(stemness_genes, rownames(cds))
cat("Available stemness genes:", paste(genes_in_data, collapse = ", "), "\n")

if (length(genes_in_data) > 0) {
  # Calculate stemness score if not already present
  if (!"Stemness_Score1" %in% colnames(colData(cds))) {
    seurat_obj <- AddModuleScore(seurat_obj, features = list(genes_in_data), name = "Stemness_Score")
    colData(cds)$Stemness_Score <- seurat_obj$Stemness_Score1
  } else {
    colData(cds)$Stemness_Score <- seurat_obj$Stemness_Score1
  }

  # Find the cluster with the highest average stemness score
  root_cluster_info <- as.data.frame(colData(cds)) %>%
    group_by(assigned_cluster) %>%
    summarise(mean_stemness = mean(Stemness_Score, na.rm = TRUE),
              n_cells = n()) %>%
    arrange(desc(mean_stemness))
  
  cat("Stemness scores by cluster:\n")
  print(root_cluster_info)
  
  root_cluster <- root_cluster_info$assigned_cluster[1]
  cat("Identified root cluster based on max stemness score:", root_cluster, "\n")

  # Get cells from the root cluster
  root_cells <- rownames(colData(cds))[colData(cds)$assigned_cluster == root_cluster]
  cat("Number of root cells:", length(root_cells), "\n")
  
} else {
  cat("WARNING: No stemness marker genes found. Using cluster 0 as root.\n")
  root_cluster <- levels(seurat_obj$seurat_clusters)[1]
  root_cells <- rownames(colData(cds))[colData(cds)$assigned_cluster == root_cluster]
}

# --- 10. Order Cells (Calculate Pseudotime) ---
cat("Step 10: Ordering cells and calculating pseudotime...\n")
tryCatch({
  cds <- order_cells(cds, root_cells = root_cells)
  colData(cds)$pseudotime <- pseudotime(cds)
  cat("Pseudotime calculated successfully.\n")
  cat("Pseudotime range:", range(pseudotime(cds), na.rm = TRUE), "\n")
}, error = function(e) {
  cat("Error calculating pseudotime:", e$message, "\n")
  cat("This might be due to disconnected graph components.\n")
  # Set pseudotime to NA for all cells if calculation fails
  colData(cds)$pseudotime <- rep(NA, ncol(cds))
})

# --- 11. Generate and Save Plots ---
cat("Step 11: Generating and saving trajectory plots...\n")
pdf(output_plot_file, width = 14, height = 10)

# Plot 1: Trajectory by Cluster
tryCatch({
  p_clusters <- plot_cells(cds, 
                          color_cells_by = "assigned_cluster", 
                          label_groups_by_cluster = TRUE,
                          label_leaves = FALSE, 
                          label_branch_points = FALSE, 
                          cell_size = 0.5,
                          alpha = 0.7) +
                ggtitle("Neftel Trajectory by Cluster") +
                theme_minimal() +
                theme(legend.position = "right")
  print(p_clusters)
}, error = function(e) {
  cat("Error plotting trajectory by cluster:", e$message, "\n")
})

# Plot 2: Trajectory by Pseudotime
if (!all(is.na(colData(cds)$pseudotime))) {
  tryCatch({
    p_pseudotime <- plot_cells(cds, 
                              color_cells_by = "pseudotime", 
                              label_cell_groups = FALSE,
                              label_leaves = FALSE, 
                              label_branch_points = FALSE, 
                              cell_size = 0.5,
                              alpha = 0.7) +
                    ggtitle("Neftel Trajectory by Pseudotime") +
                    scale_color_viridis_c(name = "Pseudotime") +
                    theme_minimal()
    print(p_pseudotime)
  }, error = function(e) {
    cat("Error plotting trajectory by pseudotime:", e$message, "\n")
  })
} else {
  cat("Skipping pseudotime plot due to calculation failure.\n")
}

# Plot 3: Stemness Score
if ("Stemness_Score" %in% colnames(colData(cds))) {
  tryCatch({
    p_stemness <- plot_cells(cds, 
                            color_cells_by = "Stemness_Score", 
                            label_cell_groups = FALSE,
                            label_leaves = FALSE, 
                            label_branch_points = FALSE, 
                            cell_size = 0.5,
                            alpha = 0.7) +
                  ggtitle("Neftel Trajectory by Stemness Score") +
                  scale_color_viridis_c(name = "Stemness\nScore") +
                  theme_minimal()
    print(p_stemness)
  }, error = function(e) {
    cat("Error plotting stemness score:", e$message, "\n")
  })
}

# Plots 4-6: UPR Pathway Scores
upr_scores <- c("IRE1_score1", "PERK_score1", "ATF6_score1")
upr_names <- c("IRE1", "PERK", "ATF6")

for (i in seq_along(upr_scores)) {
  if (upr_scores[i] %in% colnames(colData(cds))) {
    tryCatch({
      p_upr <- plot_cells(cds, 
                         color_cells_by = upr_scores[i], 
                         label_cell_groups = FALSE,
                         label_leaves = FALSE, 
                         label_branch_points = FALSE, 
                         cell_size = 0.5,
                         alpha = 0.7) +
               ggtitle(paste("Neftel Trajectory by", upr_names[i], "Score")) +
               scale_color_viridis_c(name = paste(upr_names[i], "\nScore")) +
               theme_minimal()
      print(p_upr)
    }, error = function(e) {
      cat("Error plotting", upr_names[i], "score:", e$message, "\n")
    })
  } else {
    cat("UPR score", upr_scores[i], "not found in dataset.\n")
  }
}

dev.off()
cat("Plots saved to:", output_plot_file, "\n")

# --- 12. Summary Statistics ---
cat("Step 12: Generating summary statistics...\n")
cat("\n=== TRAJECTORY ANALYSIS SUMMARY ===\n")
cat("Dataset dimensions:", dim(cds), "\n")
cat("Number of clusters:", length(unique(colData(cds)$assigned_cluster)), "\n")
cat("Root cluster:", root_cluster, "\n")
cat("Pseudotime calculated:", !all(is.na(colData(cds)$pseudotime)), "\n")

if (!all(is.na(colData(cds)$pseudotime))) {
  cat("Pseudotime statistics:\n")
  cat("  Min:", min(colData(cds)$pseudotime, na.rm = TRUE), "\n")
  cat("  Max:", max(colData(cds)$pseudotime, na.rm = TRUE), "\n")
  cat("  Mean:", mean(colData(cds)$pseudotime, na.rm = TRUE), "\n")
  cat("  Cells with valid pseudotime:", sum(!is.na(colData(cds)$pseudotime)), "\n")
}

# --- 13. Save Final Monocle CDS Object ---
cat("Step 13: Saving the final Monocle CDS object...\n")
saveRDS(cds, file = output_cds_file)
cat("Analysis complete. Final CDS object saved to '", output_cds_file, "'.\n")

# Stop logging
sink()

cat("Trajectory analysis completed successfully!\n")
cat("Check the log file for detailed information:", output_log_file, "\n")