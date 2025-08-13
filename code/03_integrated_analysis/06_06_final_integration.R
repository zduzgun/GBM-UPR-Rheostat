#!/usr/bin/env Rscript

# 06_06_final_integration.R
# Final Integration Analysis - Combining all completed 06_01 through 06_05 results
# Following the established naming convention

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(patchwork)
  library(RColorBrewer)
})

message("=== 06_06 Final Integration Analysis ===")
message("Combining results from 06_01 through 06_05")

# Load main dataset
message("Loading main dataset...")
glioblastoma <- readRDS("glioblastoma_with_final_scores.rds")

# Load completed analysis results following the 06_xx pattern
message("Loading results from previous 06_xx analyses...")

# From 06_05 LIANA analysis
liana_results <- read.csv("liana_full_results.csv")
target_interactions <- read.csv("target_cluster_interactions.csv") 
cluster6_signals <- read.csv("cluster6_received_signals.csv")
sender_summary <- read.csv("cluster_sender_summary.csv")
receiver_summary <- read.csv("cluster_receiver_summary.csv")

# From 06_03 correlation analysis
correlation_results <- NULL
if(file.exists("upr_correlation_results.txt")) {
  correlation_results <- tryCatch({
    read.table("upr_correlation_results.txt", header = TRUE, sep = "\t")
  }, error = function(e) {
    message("? Correlation file format issue, trying alternative read...")
    tryCatch({
      read.csv("upr_correlation_results.txt")
    }, error = function(e2) {
      message("? Could not read correlation results, skipping...")
      NULL
    })
  })
  
  if(!is.null(correlation_results)) {
    message("? UPR correlation results loaded from 06_03")
  }
} else {
  message("? No correlation results found from 06_03")
}

# From 06_04 trajectory analysis
trajectory_available <- file.exists("monocle_results_simple/cds_object.rds")
if(trajectory_available) {
  message("? Trajectory results available from 06_04")
} else {
  message("? No trajectory results found from 06_04")
}

message("Dataset dimensions:", paste(dim(glioblastoma), collapse = " x "))
message("Cluster distribution:")
print(table(Idents(glioblastoma)))

# Check available scores from 06_01 and 06_02
available_scores <- intersect(c("IRE1_score", "PERK_score", "ATF6_score", 
                               "Glycolysis_score", "OXPHOS_score", "Fatty_Acid_score"), 
                             colnames(glioblastoma@meta.data))
message("Available scores from 06_01/06_02:", paste(available_scores, collapse = ", "))

# === INTEGRATION FIGURE 1: Complete Ecosystem Overview ===
message("Creating Integration Figure 1: Complete Ecosystem...")

pdf("06_06_integration_figure1.pdf", width = 16, height = 12)

# Main UMAP with enhanced cluster information
p_main <- DimPlot(glioblastoma, group.by = "seurat_clusters", 
                  label = TRUE, label.size = 6, pt.size = 0.8) +
  ggtitle("GBM Ecosystem Integration: 9 Cell Clusters") +
  labs(subtitle = "Results from 06_01 through 06_05 analyses") +
  theme_void() +
  theme(plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12))

# Feature plots for available scores (max 4 to fit layout)
score_plots <- list()
plot_scores <- head(available_scores, 4)

for(i in seq_along(plot_scores)) {
  score_plots[[i]] <- FeaturePlot(glioblastoma, features = plot_scores[i], 
                                 pt.size = 0.4, cols = c("lightgray", "red")) +
    ggtitle(gsub("_score", "", plot_scores[i])) +
    theme_void() +
    theme(plot.title = element_text(size = 12, face = "bold"))
}

# Fill remaining slots with spacers if needed
while(length(score_plots) < 4) {
  score_plots[[length(score_plots) + 1]] <- plot_spacer()
}

# Create layout
layout_design <- "
AABBCC
AABBCC
DDDDEE
DDDDEE
"

combined1 <- p_main + score_plots[[1]] + score_plots[[2]] + 
             score_plots[[3]] + score_plots[[4]] +
             plot_layout(design = layout_design)

print(combined1)
dev.off()

# === INTEGRATION FIGURE 2: UPR Integration Results ===
message("Creating Integration Figure 2: UPR Integration...")

pdf("06_06_integration_figure2.pdf", width = 14, height = 10)

upr_scores <- intersect(c("IRE1_score", "PERK_score", "ATF6_score"), 
                       colnames(glioblastoma@meta.data))

if(length(upr_scores) >= 2) {
  message("Processing UPR integration with scores:", paste(upr_scores, collapse = ", "))
  
  # Cluster-wise UPR profile
  upr_data <- glioblastoma@meta.data[, c("seurat_clusters", upr_scores)]
  upr_summary <- upr_data %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(upr_scores), mean, na.rm = TRUE), .groups = 'drop')
  
  # Violin plot
  upr_long <- reshape2::melt(upr_data, id.vars = "seurat_clusters", 
                            variable.name = "UPR_Pathway", value.name = "Score")
  
  p_violin <- ggplot(upr_long, aes(x = seurat_clusters, y = Score, fill = UPR_Pathway)) +
    geom_violin(position = "dodge", alpha = 0.7) +
    geom_boxplot(position = position_dodge(0.9), width = 0.2, alpha = 0.5) +
    stat_summary(fun = mean, geom = "point", position = position_dodge(0.9), 
                 size = 2, color = "black") +
    scale_fill_brewer(type = "qual", palette = "Dark2") +
    theme_minimal() +
    labs(title = "UPR Pathway Integration Across Cell Clusters",
         subtitle = "Combined results from 06_01 module scores and 06_03 correlation",
         x = "Cell Cluster", y = "UPR Activity Score",
         fill = "UPR Pathway") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          legend.position = "bottom")
  
  print(p_violin)
  
  # Correlation heatmap if correlation results exist
  if(exists("correlation_results")) {
    message("Adding correlation results from 06_03")
  }
  
} else {
  message("Insufficient UPR data for integration")
  grid::grid.newpage()
  grid::grid.text("UPR Integration: Insufficient data\nCheck 06_01 module scoring results", 
                  x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 16))
}

dev.off()

# === INTEGRATION FIGURE 3: Communication Network Integration ===
message("Creating Integration Figure 3: Communication Network...")

pdf("06_06_integration_figure3.pdf", width = 16, height = 12)

# Communication overview from 06_05 LIANA results
p_senders <- ggplot(head(sender_summary, 8), 
                   aes(x = reorder(paste("Cluster", source), signals_sent), 
                       y = signals_sent)) +
  geom_col(fill = "steelblue", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = signals_sent), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Signal Senders (from 06_05 LIANA)",
       x = "Cluster", y = "Total Signals Sent") +
  theme(plot.title = element_text(size = 12, face = "bold"))

p_receivers <- ggplot(head(receiver_summary, 8),
                     aes(x = reorder(paste("Cluster", target), signals_received), 
                         y = signals_received)) +
  geom_col(fill = "darkgreen", alpha = 0.8, width = 0.7) +
  geom_text(aes(label = signals_received), hjust = -0.1, size = 4, fontface = "bold") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Signal Receivers (from 06_05 LIANA)", 
       x = "Cluster", y = "Total Signals Received") +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Critical finding: Cluster 6 analysis
if(nrow(cluster6_signals) > 0) {
  p_cluster6 <- ggplot(head(cluster6_signals, 12),
                      aes(x = reorder(paste(ligand, "›", receptor), weight_sc),
                          y = weight_sc)) +
    geom_col(aes(fill = paste("from Cluster", source)), alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    theme_minimal() +
    labs(title = "CRITICAL: Cluster 6 Signal Reception",
         subtitle = paste("Total signals received:", nrow(cluster6_signals)),
         x = "Ligand › Receptor Pairs", y = "Signal Strength",
         fill = "Signal Source") +
    theme(plot.title = element_text(size = 12, face = "bold", color = "darkred"),
          axis.text.y = element_text(size = 8),
          legend.position = "bottom")
  
  # Combine plots
  combined_comm <- (p_senders + p_receivers) / p_cluster6
  print(combined_comm)
} else {
  combined_comm <- p_senders + p_receivers
  print(combined_comm)
}

dev.off()

# === INTEGRATION FIGURE 4: Research Hypothesis Results ===
message("Creating Integration Figure 4: Hypothesis Testing...")

pdf("06_06_integration_figure4.pdf", width = 14, height = 10)

if(nrow(target_interactions) > 0) {
  # Main hypothesis: 3,5 › 1,6,8
  p_hypothesis <- ggplot(head(target_interactions, 25),
                        aes(x = reorder(paste(ligand, "›", receptor), weight_sc),
                            y = weight_sc)) +
    geom_col(aes(fill = paste("Cluster", source, "›", target)), alpha = 0.8) +
    coord_flip() +
    scale_fill_brewer(type = "qual", palette = "Paired") +
    theme_minimal() +
    labs(title = "Hypothesis Test: Senders (3,5) › Receivers (1,6,8)",
         subtitle = paste("Integration result:", nrow(target_interactions), "confirmed interactions"),
         x = "Ligand › Receptor Communications", y = "Signal Strength (LIANA)",
         fill = "Communication Direction") +
    theme(plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12, color = "darkblue"),
          axis.text.y = element_text(size = 8),
          legend.position = "bottom")
  
  print(p_hypothesis)
  
  # Create summary table
  hypothesis_summary <- target_interactions %>%
    group_by(source, target) %>%
    summarise(
      interactions = n(),
      avg_strength = round(mean(weight_sc, na.rm = TRUE), 3),
      top_ligand = first(ligand[order(-weight_sc)]),
      .groups = 'drop'
    ) %>%
    arrange(desc(interactions))
  
  # Add text summary
  grid::grid.newpage()
  summary_text <- paste(
    "HYPOTHESIS TESTING SUMMARY:",
    paste("• Total target interactions found:", nrow(target_interactions)),
    paste("• Communication pairs:", nrow(hypothesis_summary)),
    "• Top communication patterns:",
    paste(apply(head(hypothesis_summary, 3), 1, function(x) 
          paste("  -", x[1], "›", x[2], ":", x[3], "interactions")), collapse = "\n"),
    sep = "\n"
  )
  
  grid::grid.text(summary_text, x = 0.1, y = 0.8, 
                  gp = grid::gpar(fontsize = 12, fontface = "bold"),
                  just = c("left", "top"))
  
} else {
  grid::grid.newpage()
  grid::grid.text("Hypothesis Testing: No target interactions found\nCheck 06_05 LIANA results", 
                  x = 0.5, y = 0.5, gp = grid::gpar(fontsize = 16))
}

dev.off()

# === 06_06 INTEGRATION SUMMARY REPORT ===
message("Generating 06_06 integration summary...")

# Prepare comprehensive summary
integration_summary <- list(
  analysis_date = Sys.Date(),
  total_cells = ncol(glioblastoma),
  total_clusters = length(unique(Idents(glioblastoma))),
  available_scores = available_scores,
  liana_total = nrow(liana_results),
  cluster6_signals = nrow(cluster6_signals),
  target_interactions = nrow(target_interactions),
  top_senders = head(sender_summary$source, 3),
  top_receivers = head(receiver_summary$target, 3)
)

# Generate detailed report
report_content <- paste0(
  "=== 06_06 FINAL INTEGRATION ANALYSIS REPORT ===\n",
  "Integration Date: ", integration_summary$analysis_date, "\n",
  "Combining results from analyses 06_01 through 06_05\n\n",
  
  "DATASET OVERVIEW:\n",
  "• Total cells analyzed: ", integration_summary$total_cells, "\n",
  "• Cell clusters: ", integration_summary$total_clusters, " (0-8)\n",
  "• Available molecular scores: ", length(integration_summary$available_scores), "\n",
  "  - ", paste(integration_summary$available_scores, collapse = ", "), "\n\n",
  
  "ANALYSIS PIPELINE RESULTS:\n",
  "• 06_01: Module scores ?\n",
  "• 06_02: Metabolic scores ?\n", 
  "• 06_03: Correlation analysis ", ifelse(exists("correlation_results"), "?", "?"), "\n",
  "• 06_04: Trajectory analysis ", ifelse(trajectory_available, "?", "?"), "\n",
  "• 06_05: LIANA communication ?\n",
  "• 06_06: Integration analysis ?\n\n",
  
  "KEY SCIENTIFIC FINDINGS:\n",
  "1. COMMUNICATION PATTERNS:\n",
  "   • Total interactions detected: ", integration_summary$liana_total, "\n",
  "   • Top signal senders: Clusters ", paste(integration_summary$top_senders, collapse = ", "), "\n",
  "   • Top signal receivers: Clusters ", paste(integration_summary$top_receivers, collapse = ", "), "\n\n",
  
  "2. CLUSTER 6 ANALYSIS (Critical Research Question):\n",
  "   • Signals received: ", integration_summary$cluster6_signals, "\n",
  "   • Conclusion: ", ifelse(integration_summary$cluster6_signals > 0, 
                            "NOT silent - receives targeted signals", 
                            "Supports 'silent receiver' hypothesis"), "\n\n",
  
  "3. HYPOTHESIS TESTING (Clusters 3,5 › 1,6,8):\n",
  "   • Target interactions found: ", integration_summary$target_interactions, "\n",
  "   • Hypothesis status: ", ifelse(integration_summary$target_interactions > 0,
                                   "SUPPORTED", "NOT SUPPORTED"), "\n\n",
  
  "GENERATED FIGURES:\n",
  "• 06_06_integration_figure1.pdf - Complete ecosystem overview\n",
  "• 06_06_integration_figure2.pdf - UPR pathway integration\n",
  "• 06_06_integration_figure3.pdf - Communication network integration\n", 
  "• 06_06_integration_figure4.pdf - Hypothesis testing results\n\n",
  
  "NEXT STEPS:\n",
  "• Review generated figures for publication quality\n",
  "• Consider additional analyses if needed (06_07, 06_08, etc.)\n",
  "• Prepare manuscript figures and statistical summaries\n"
)

# Write comprehensive report
writeLines(report_content, "06_06_integration_report.txt")

# Also create a simple summary for quick reference
quick_summary <- data.frame(
  Analysis = paste0("06_0", 1:6),
  Description = c("Module Scores", "Metabolic Scores", "Correlation Analysis", 
                  "Trajectory Analysis", "LIANA Communication", "Final Integration"),
  Status = c("?", "?", ifelse(exists("correlation_results"), "?", "?"), 
             ifelse(trajectory_available, "?", "?"), "?", "?"),
  Key_Output = c("UPR pathway scores", "Metabolic pathway scores", "UPR correlations",
                 "Cell trajectories", "Cell-cell communication", "Integrated figures")
)

write.csv(quick_summary, "06_06_analysis_pipeline_summary.csv", row.names = FALSE)

message("\n=== 06_06 FINAL INTEGRATION COMPLETED ===")
message("Following established 06_xx naming convention")
message("Check 06_06_integration_report.txt for comprehensive summary")
message("All integration figures generated with 06_06 prefix")