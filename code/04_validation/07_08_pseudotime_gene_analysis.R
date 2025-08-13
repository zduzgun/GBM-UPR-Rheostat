# ==============================================================================
# 07_08_pseudotime_gene_analysis.R (v3 - find_gene_modules Fix)
#
# Purpose: Analyze gene expression dynamics along the malignant cell trajectory.
#          FIX: Runs preprocess_cds before find_gene_modules to build PCA model.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
input_cds_file <- "07_07_malignant_monocle_cds.rds"
output_gene_table <- "results/tables/07_08_pseudotime_dependent_genes.csv"
output_plot_file <- "results/figures/07_08_pseudotime_gene_plots.pdf"
output_final_cds <- "07_08_final_cds_with_modules.rds"

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

# --- 3. Load Monocle CDS Object ---
cat("Step 3: Loading the trajectory CDS object from 07_07...\n")
cds <- readRDS(input_cds_file)
cat("CDS object loaded successfully.\n")

# --- 4. Find Genes that Vary Along Pseudotime ---
cat("Step 4: Identifying genes with expression dependent on pseudotime...\n")
graph_results <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

significant_genes <- graph_results %>%
  arrange(q_value) %>%
  filter(q_value < 0.05)

cat("  - Found", nrow(significant_genes), "genes that vary significantly along the trajectory (q < 0.05).\n")
write.csv(significant_genes, file = output_gene_table, row.names = FALSE)

# --- 5. Visualize Top Pseudotime-Dependent Genes ---
cat("Step 5: Visualizing expression of top 6 significant genes...\n")
top_6_genes <- head(significant_genes, 6)

pdf(output_plot_file, width = 12, height = 8)

p_top_genes_umap <- plot_cells(cds,
                               genes = top_6_genes$gene_short_name,
                               show_trajectory_graph = TRUE,
                               label_cell_groups = FALSE,
                               label_leaves = FALSE,
                               cell_size = 0.5) +
  ggtitle("Expression of Top 6 Pseudotime-Dependent Genes")
print(p_top_genes_umap)

p_top_genes_pseudotime <- plot_genes_in_pseudotime(cds[top_6_genes$gene_short_name, ],
                                                   color_cells_by = "pseudotime",
                                                   min_expr = 0.5) +
  ggtitle("Expression vs. Pseudotime for Top 6 Genes")
print(p_top_genes_pseudotime)

# --- 6. Re-visit UPR & Reostat Scores Along the Malignant Trajectory ---
cat("Step 6: Analyzing UPR pathway scores along the trajectory...\n")
upr_scores_to_plot <- c("IRE1_score", "PERK_score", "ATF6_score", "UPR_Reostat_SD")
available_scores <- upr_scores_to_plot[upr_scores_to_plot %in% colnames(colData(cds))]

if (length(available_scores) > 0) {
  cat("  - Plotting available UPR scores:", paste(available_scores, collapse = ", "), "\n")
  p_upr_trajectory <- plot_cells(cds,
                                 color_cells_by = available_scores[1],
                                 show_trajectory_graph = TRUE,
                                 label_cell_groups = FALSE,
                                 cell_size = 0.5) +
    ggtitle(paste(available_scores[1], "Along Malignant Trajectory")) +
    scale_color_viridis_c()
  print(p_upr_trajectory)
  
  cds$pseudotime <- pseudotime(cds)
  p_upr_over_time <- ggplot(as.data.frame(colData(cds)), aes(x = pseudotime, y = .data[[available_scores[1]]])) +
    geom_point(alpha = 0.2, size=0.8) +
    geom_smooth(method = "loess", color = "red") +
    ggtitle(paste(available_scores[1], "vs. Pseudotime")) +
    theme_minimal()
  print(p_upr_over_time)

} else {
  cat("  - WARNING: No UPR score columns found in the CDS object metadata.\n")
}

# --- 7. Group Significant Genes into Co-expression Modules ---
cat("Step 7: Clustering pseudotime-dependent genes into modules...\n")

# ================================ DÜZELTME BURADA =================================
# Çözüm: find_gene_modules'dan önce preprocess_cds'i çalýþtýrarak
#        fonksiyonun ihtiyaç duyduðu PCA modelini cds nesnesi içinde oluþtur.
cat("  - Pre-processing CDS with PCA to build the required model for find_gene_modules...\n")
cds <- preprocess_cds(cds, method = "PCA")
# ===============================================================================

# Bu gruplama, ifade desenleri benzer olan genleri bir araya getirir.
gene_modules <- find_gene_modules(cds[significant_genes$gene_short_name, ],
                                  resolution = 0.01,
                                  cores = 4)

# Modül bilgisini çizim için CDS nesnesine ekle
colData(cds)$gene_module <- gene_modules[match(rownames(colData(cds)), names(gene_modules))]
write.csv(gene_modules, file = "results/tables/07_08_gene_modules.csv")

p_gene_modules <- plot_cells(cds,
                             color_cells_by = "gene_module",
                             show_trajectory_graph = TRUE,
                             label_cell_groups = FALSE,
                             label_leaves = FALSE,
                             cell_size = 0.5) +
  ggtitle("Gene Modules Along Malignant Trajectory")
print(p_gene_modules)

dev.off()

# --- 8. Save Final Object ---
cat("Step 8: Saving the final CDS object with module information...\n")
saveRDS(cds, file = output_final_cds)
cat("Analysis complete. Final objects and plots saved successfully.\n")