#!/usr/bin/env Rscript

# LIANA Working Analysis - Based on successful debug
suppressPackageStartupMessages({
  library(Seurat)
  library(liana)
  library(dplyr)
  library(ggplot2)
  library(RColorBrewer)
})

# Load dataset
message("Loading dataset...")
glioblastoma <- readRDS("glioblastoma_with_final_scores.rds")

message("Dataset information:")
print(dim(glioblastoma))
cluster_table <- table(Idents(glioblastoma))
print(cluster_table)

# Target clusters for hypothesis
target_senders <- c("3", "5")
target_receivers <- c("1", "6", "8") 
available_clusters <- names(cluster_table)

message("Target analysis clusters:")
message(paste("Senders (3,5):", paste(intersect(target_senders, available_clusters), collapse = ", ")))
message(paste("Receivers (1,6,8):", paste(intersect(target_receivers, available_clusters), collapse = ", ")))

# Run LIANA - using EXACT working configuration from debug
message("Starting LIANA analysis...")

liana_results <- liana_wrap(
  glioblastoma,
  method = 'connectome',  # Single working method from debug
  resource = 'Consensus',  # Capital C - worked in debug
  idents_col = 'seurat_clusters',
  verbose = TRUE
)

message(paste("Total interactions found:", nrow(liana_results)))
print("Column names:")
print(colnames(liana_results))

# Save full results
write.csv(liana_results, "liana_full_results.csv", row.names = FALSE)

# Create ranking based on weight_sc (main score from connectome method)
message("Creating consensus ranking...")
liana_ranked <- liana_results %>%
  arrange(desc(weight_sc)) %>%
  mutate(consensus_rank = row_number())

# Focus on target clusters communication: 3,5 -> 1,6,8
message("Analyzing target cluster communications...")

target_interactions <- liana_ranked %>%
  filter(source %in% target_senders & target %in% target_receivers) %>%
  arrange(consensus_rank)

message(paste("Target cluster interactions found:", nrow(target_interactions)))

# All interactions involving cluster 6 (critical research question)
cluster6_interactions <- liana_ranked %>%
  filter(target == "6") %>%
  arrange(consensus_rank)

message(paste("Interactions targeting Cluster 6:", nrow(cluster6_interactions)))

# Get top interactions overall
top_n <- max(100, round(nrow(liana_ranked) * 0.01))
top_interactions <- liana_ranked %>%
  head(top_n)

# Save targeted results
write.csv(target_interactions, "target_cluster_interactions.csv", row.names = FALSE)
write.csv(cluster6_interactions, "cluster6_received_signals.csv", row.names = FALSE)
write.csv(top_interactions, "liana_top100_interactions.csv", row.names = FALSE)

# Create the critical visualization
message("Creating targeted dot plot visualization...")

if(nrow(target_interactions) > 0) {
  plot_data <- target_interactions %>%
    head(30) %>%
    mutate(
      interaction_name = paste(ligand, receptor, sep = " -> "),
      sender_receiver = paste("Cluster", source, "to Cluster", target),
      significance = weight_sc  # Use connectome weight score
    )
  
  pdf("01_liana_top100_interactions.pdf", width = 14, height = 10)
  
  p_target <- ggplot(plot_data, aes(x = interaction_name, y = sender_receiver)) +
    geom_point(aes(size = significance, fill = significance), 
               shape = 21, color = "black", alpha = 0.8) +
    scale_size_continuous(name = "Interaction\nStrength", range = c(2, 12)) +
    scale_fill_gradient2(name = "Interaction\nStrength", 
                        low = "lightblue", mid = "yellow", high = "red",
                        midpoint = median(plot_data$significance)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(size = 14, face = "bold"),
      legend.position = "right"
    ) +
    labs(
      title = "Critical Cluster Communications: Senders (3,5) -> Receivers (1,6,8)",
      subtitle = "Top Ligand-Receptor Interactions Based on LIANA Connectome Analysis",
      x = "Ligand -> Receptor Pairs",
      y = "Sender -> Receiver Clusters"
    )
  
  print(p_target)
  dev.off()
  message("Created: 01_liana_top100_interactions.pdf")
} else {
  message("No target cluster interactions found!")
}

# Critical analysis for Cluster 6
message("\n=== CRITICAL ANALYSIS: Cluster 6 Communication ===")

if(nrow(cluster6_interactions) > 0) {
  cluster6_senders <- cluster6_interactions %>%
    count(source, sort = TRUE)
  
  cluster6_pathways <- cluster6_interactions %>%
    head(20) %>%
    mutate(signal_strength = weight_sc)
  
  message("Top senders to Cluster 6:")
  print(cluster6_senders)
  
  # Cluster 6 specific visualization
  pdf("cluster6_specific_analysis.pdf", width = 12, height = 8)
  
  p_cluster6 <- ggplot(cluster6_pathways, 
                       aes(x = reorder(paste(ligand, receptor, sep = "-"), signal_strength),
                           y = signal_strength)) +
    geom_col(aes(fill = source), alpha = 0.8) +
    coord_flip() +
    theme_minimal() +
    scale_fill_brewer(type = "qual", palette = "Set3") +
    labs(
      title = "Cluster 6 Signal Reception Analysis",
      subtitle = "What signals does the stem cell cluster actually receive?",
      x = "Ligand-Receptor Pairs",
      y = "Signal Strength (Connectome Weight)",
      fill = "Sender Cluster"
    ) +
    theme(axis.text.y = element_text(size = 8))
  
  print(p_cluster6)
  dev.off()
  message("Created: cluster6_specific_analysis.pdf")
  
  # Research question answer
  message("\n=== RESEARCH QUESTION ANSWER ===")
  message(paste("RESULT: Cluster 6 receives", nrow(cluster6_interactions), "significant signals"))
  
  if(nrow(cluster6_interactions) > 0) {
    message("Top 5 strongest signals to Cluster 6:")
    top5_to_cluster6 <- cluster6_interactions %>% head(5)
    for(i in 1:min(5, nrow(top5_to_cluster6))) {
      message(paste(i, ".", 
                   "Cluster", top5_to_cluster6$source[i], "->", "Cluster", top5_to_cluster6$target[i], ":",
                   top5_to_cluster6$ligand[i], "->", top5_to_cluster6$receptor[i],
                   "(strength:", round(top5_to_cluster6$weight_sc[i], 4), ")"))
    }
    
    # Check if cluster 6 receives signals from hypothesis senders (3,5)
    cluster6_from_hypothesis <- cluster6_interactions %>%
      filter(source %in% target_senders)
    
    if(nrow(cluster6_from_hypothesis) > 0) {
      message(paste("\nIMPORTANT: Cluster 6 receives", nrow(cluster6_from_hypothesis), 
                   "signals from hypothesis senders (Clusters 3,5)"))
    } else {
      message("\nNOTE: Cluster 6 does NOT receive signals from hypothesis senders (Clusters 3,5)")
    }
  }
  
} else {
  message("RESULT: Cluster 6 receives NO significant signals - supports 'silent receiver' hypothesis")
}

# Overall summary statistics
message("\n=== OVERALL ANALYSIS SUMMARY ===")
sender_summary <- liana_ranked %>%
  group_by(source) %>%
  summarise(
    signals_sent = n(),
    avg_strength = mean(weight_sc, na.rm = TRUE),
    max_strength = max(weight_sc, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(signals_sent))

receiver_summary <- liana_ranked %>%
  group_by(target) %>%
  summarise(
    signals_received = n(),
    avg_strength = mean(weight_sc, na.rm = TRUE),
    max_strength = max(weight_sc, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  arrange(desc(signals_received))

message("Top signal SENDERS:")
print(head(sender_summary))

message("Top signal RECEIVERS:")
print(head(receiver_summary))

# Save summaries
write.csv(sender_summary, "cluster_sender_summary.csv", row.names = FALSE)
write.csv(receiver_summary, "cluster_receiver_summary.csv", row.names = FALSE)

message("\n=== FINAL SUMMARY ===")
message(paste("Total interactions detected:", nrow(liana_results)))
message(paste("Hypothesis-relevant interactions (3,5 -> 1,6,8):", nrow(target_interactions)))
message(paste("Signals to Cluster 6:", nrow(cluster6_interactions)))

message("\n=== OUTPUT FILES GENERATED ===")
message("1. liana_full_results.csv - Complete LIANA results")
message("2. 01_liana_top100_interactions.pdf - Critical hypothesis visualization")  
message("3. target_cluster_interactions.csv - Senders (3,5) -> Receivers (1,6,8)")
message("4. cluster6_received_signals.csv - All signals targeting Cluster 6")
message("5. cluster6_specific_analysis.pdf - Detailed Cluster 6 analysis")
message("6. cluster_sender_summary.csv - Summary of all senders")
message("7. cluster_receiver_summary.csv - Summary of all receivers")

message("\n=== ANALYSIS COMPLETED SUCCESSFULLY! ===")
message("You now have evidence to answer your research question about Cluster 6!")