# ==============================================================================
# 06_04_trajectory_SIMPLE_ROBUST.R
#
# Purpose: Basit ama güvenilir trajectory analizi - hata-free
# Strategy: Partition problemini bypass ederek doðrudan trajectory
# ==============================================================================

# --- 1. Load Libraries ---
suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(viridis)
})

cat("Starting Simple but Robust Trajectory Analysis...\n")

# --- 2. Load Data ---
input_rds_file <- "glioblastoma_with_final_scores.rds"
seurat_data <- readRDS(input_rds_file)

if (is.list(seurat_data) && !is(seurat_data, "Seurat")) {
  seurat_obj_index <- which(sapply(seurat_data, function(x) is(x, "Seurat")))
  seurat_obj <- seurat_data[[seurat_obj_index[1]]]
} else {
  seurat_obj <- seurat_data
}

cat("Loaded Seurat object with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n")

# --- 3. Create CDS ---
expression_matrix <- GetAssayData(seurat_obj, assay = 'RNA', layer = 'counts')
cell_metadata <- seurat_obj@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(expression_matrix), 
                           row.names = rownames(expression_matrix))
cds <- new_cell_data_set(expression_matrix, cell_metadata = cell_metadata, gene_metadata = gene_metadata)

# --- 4. Preprocessing ---
cat("Preprocessing CDS object...\n")
cds <- preprocess_cds(cds, num_dim = 50)
reducedDims(cds)[['UMAP']] <- seurat_obj@reductions$umap@cell.embeddings
colData(cds)$assigned_cluster <- seurat_obj$seurat_clusters

# --- 5. Simple Clustering and Graph Learning ---
cat("Learning graph structure...\n")
cds <- cluster_cells(cds, reduction_method = "UMAP")

# CRITICAL FIX: Learn graph WITHOUT using partitions
cds <- learn_graph(cds, use_partition = FALSE)  # This bypasses partition problems

# --- 6. Add Stemness Scores ---
stemness_genes <- c('SOX2', 'NES', 'PROM1', 'OLIG2', 'MYCN', 'MAGEA3')
final_stemness_genes <- intersect(stemness_genes, rownames(cds))

cat("Using stemness genes:", paste(final_stemness_genes, collapse = ", "), "\n")

if (!"Stemness_Score1" %in% colnames(seurat_obj@meta.data)) {
    seurat_obj <- AddModuleScore(seurat_obj, features = list(final_stemness_genes), name = "Stemness_Score")
}
colData(cds)$Stemness_Score1 <- seurat_obj@meta.data$Stemness_Score1

# Transfer UPR scores if available
upr_scores <- c("IRE1_Score1", "PERK_Score1", "ATF6_Score1")
for (score in upr_scores) {
  if (score %in% colnames(seurat_obj@meta.data)) {
    colData(cds)[[score]] <- seurat_obj@meta.data[[score]]
    cat("Transferred", score, "\n")
  }
}

# --- 7. Simple Root Selection ---
cat("Selecting root cluster based on stemness...\n")
cluster_stemness <- as.data.frame(colData(cds)) %>%
  group_by(assigned_cluster) %>%
  summarise(
    cell_count = n(),
    mean_stemness = mean(Stemness_Score1, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  filter(cell_count >= 15) %>%  # At least 15 cells
  arrange(desc(mean_stemness))

cat("Cluster analysis:\n")
print(cluster_stemness)

if (nrow(cluster_stemness) == 0) {
  stop("No clusters with sufficient cells found!")
}

root_cluster <- cluster_stemness$assigned_cluster[1]
cat("Selected root cluster:", as.character(root_cluster), "\n")

# --- 8. Order Cells ---
cat("Ordering cells along trajectory...\n")
root_cells <- colnames(cds)[colData(cds)$assigned_cluster == root_cluster]

if (length(root_cells) == 0) {
  stop("No root cells found!")
}

cat("Found", length(root_cells), "root cells\n")

# Order cells
tryCatch({
  cds <- order_cells(cds, root_cells = root_cells)
  
  # Check results
  pt_values <- pseudotime(cds)
  valid_pt <- sum(!is.na(pt_values))
  
  cat("SUCCESS: Pseudotime calculated for", valid_pt, "/", ncol(cds), "cells\n")
  cat("Pseudotime range:", round(min(pt_values, na.rm = TRUE), 2), "-", 
      round(max(pt_values, na.rm = TRUE), 2), "\n")
  
}, error = function(e) {
  cat("ERROR in ordering cells:", e$message, "\n")
  stop("Trajectory calculation failed")
})

# --- 9. Generate Comprehensive Plots ---
cat("Generating comprehensive visualizations...\n")
output_file <- "glioblastoma_trajectory_SIMPLE_ROBUST.pdf"

pdf(output_file, width = 14, height = 10)

# Plot 1: Clusters
p1 <- plot_cells(cds, 
                color_cells_by = "assigned_cluster", 
                label_groups_by_cluster = TRUE, 
                cell_size = 0.6,
                show_trajectory_graph = TRUE) +
  ggtitle("GBM Trajectory: Seurat Clusters") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p1)

# Plot 2: Pseudotime
p2 <- plot_cells(cds, 
                color_cells_by = "pseudotime", 
                cell_size = 0.6,
                show_trajectory_graph = TRUE) +
  ggtitle("GBM Trajectory: Pseudotime") +
  scale_color_viridis_c(name = "Pseudotime", na.value = "grey90") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p2)

# Plot 3: Stemness Score
p3 <- plot_cells(cds, 
                color_cells_by = "Stemness_Score1", 
                cell_size = 0.6,
                show_trajectory_graph = TRUE) +
  ggtitle("GBM Trajectory: Stemness Score") +
  scale_color_viridis_c(name = "Stemness\nScore") +
  theme_minimal() +
  theme(legend.position = "bottom")
print(p3)

# Plot 4-6: UPR Pathways
for (score in upr_scores) {
  if (score %in% colnames(colData(cds))) {
    pathway_name <- gsub("_Score1", "", score)
    p_pathway <- plot_cells(cds, 
                           color_cells_by = score, 
                           cell_size = 0.6,
                           show_trajectory_graph = TRUE) +
      ggtitle(paste("GBM Trajectory:", pathway_name, "Activity")) +
      scale_color_viridis_c(name = paste(pathway_name, "\nScore")) +
      theme_minimal() +
      theme(legend.position = "bottom")
    print(p_pathway)
  }
}

# Plot 7: Stemness genes expression
if (length(final_stemness_genes) > 0) {
  p_genes <- plot_cells(cds, 
                       genes = final_stemness_genes, 
                       cell_size = 0.4,
                       show_trajectory_graph = FALSE) +
    ggtitle("Stemness Genes Expression") +
    theme_minimal()
  print(p_genes)
}

dev.off()

cat("Plots saved to:", output_file, "\n")

# --- 10. Save Results ---
cat("Saving analysis results...\n")

# Save CDS object
output_dir <- "monocle_results_simple"
if (!dir.exists(output_dir)) dir.create(output_dir)

tryCatch({
  save_monocle_objects(cds, directory = output_dir)
  cat("CDS object saved to:", output_dir, "\n")
}, error = function(e) {
  cat("WARNING: Could not save with save_monocle_objects. Saving as RDS...\n")
  saveRDS(cds, file = file.path(output_dir, "cds_with_trajectory.rds"))
})

# --- 11. Final Summary ---
cat("\n=== FINAL RESULTS SUMMARY ===\n")
cat("Total cells analyzed:", ncol(cds), "\n")
cat("Root cluster:", as.character(root_cluster), "\n")
cat("Cells with pseudotime:", sum(!is.na(pseudotime(cds))), "\n")
cat("Stemness genes used:", paste(final_stemness_genes, collapse = ", "), "\n")

# Pseudotime statistics
pt_stats <- summary(pseudotime(cds))
cat("Pseudotime statistics:\n")
print(pt_stats)

# Cluster-wise pseudotime summary
cat("\nPseudotime by cluster:\n")
cluster_pt_summary <- as.data.frame(colData(cds)) %>%
  group_by(assigned_cluster) %>%
  summarise(
    cells = n(),
    mean_pseudotime = mean(pseudotime, na.rm = TRUE),
    median_pseudotime = median(pseudotime, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(mean_pseudotime)

print(cluster_pt_summary)

cat("\n? SIMPLE TRAJECTORY ANALYSIS COMPLETED SUCCESSFULLY! ?\n")
cat("?? Results saved in:", output_dir, "\n")
cat("?? Plots saved as:", output_file, "\n")