#!/usr/bin/env Rscript
# Figure 2: UPR Rheostat Hypothesis - ASCII safe version
# Enhanced version with improved visuals and publication-quality styling

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(scales)
  library(viridis)
})

# Input file
input_rds <- "glioblastoma_with_corrected_rheostat.rds"
output_name <- "Figure2_Enhanced_UPR_Rheostat_Hypothesis"

# Helper functions
choose_first_available <- function(cn, patterns) {
  for (pt in patterns) {
    idx <- which(grepl(pt, cn, ignore.case = TRUE))
    if (length(idx) > 0) return(cn[idx[1]])
  }
  return(NA_character_)
}

get_metadata_df <- function(x) {
  if (inherits(x, "Seurat")) return(x@meta.data %>% tibble::rownames_to_column(var = "cell_id"))
  if (inherits(x, "SingleCellExperiment")) return(as.data.frame(SummarizedExperiment::colData(x)) %>% tibble::rownames_to_column(var = "cell_id"))
  stop("Object type not supported. Expected Seurat or SingleCellExperiment.")
}

# Statistical functions
calculate_correlation_with_p <- function(x, y, method = "spearman") {
  if (length(x) < 3 || length(y) < 3) return(list(cor = NA, p = NA))
  test_result <- suppressWarnings(cor.test(x, y, method = method))
  list(cor = test_result$estimate, p = test_result$p.value)
}

format_p_value <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**", 
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

# UPR rheostat state classification
classify_rheostat_state <- function(df, perk, ire1, atf6, q = 0.75, low_q = 0.25) {
  hi_perk <- df[[perk]] > quantile(df[[perk]], probs = q, na.rm = TRUE)
  hi_ire1 <- df[[ire1]] > quantile(df[[ire1]], probs = q, na.rm = TRUE)
  hi_atf6 <- df[[atf6]] > quantile(df[[atf6]], probs = q, na.rm = TRUE)
  lo_all <- df[[perk]] < quantile(df[[perk]], probs = low_q, na.rm = TRUE) &
            df[[ire1]] < quantile(df[[ire1]], probs = low_q, na.rm = TRUE) &
            df[[atf6]] < quantile(df[[atf6]], probs = low_q, na.rm = TRUE)

  state <- rep("Intermediate/Mixed", nrow(df))
  state[lo_all] <- "All-low"
  state[hi_perk & !hi_ire1 & !hi_atf6] <- "PERK-high"
  state[hi_ire1 & !hi_perk & !hi_atf6] <- "IRE1-high"
  state[hi_atf6 & !hi_perk & !hi_ire1] <- "ATF6-high"
  state[(hi_perk + hi_ire1 + hi_atf6) >= 2] <- "Multi-arm-high"
  factor(state, levels = c("All-low", "PERK-high", "IRE1-high", "ATF6-high", "Multi-arm-high", "Intermediate/Mixed"))
}

# Enhanced publication theme
theme_publication_enhanced <- function(base_size = 12) {
  theme_minimal(base_size = base_size) %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey20", fill = NA, linewidth = 0.8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
      text = element_text(colour = "grey15"),
      axis.text = element_text(colour = "grey25", size = rel(0.9)),
      axis.title = element_text(colour = "grey15", size = rel(1.0), face = "bold"),
      plot.title = element_text(colour = "grey10", size = rel(1.3), face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(colour = "grey40", size = rel(0.9), margin = margin(b = 15)),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.size = unit(0.5, "cm"),
      legend.spacing = unit(0.3, "cm"),
      legend.title = element_text(colour = "grey15", size = rel(0.9), face = "bold"),
      legend.text = element_text(colour = "grey25", size = rel(0.8)),
      legend.margin = margin(l = 10),
      plot.margin = unit(c(1.2, 1.2, 1.2, 1.2), "lines"),
      axis.line = element_line(colour = "grey30", linewidth = 0.6),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm")
    )
}

# Load and process data
cat("Loading data from:", input_rds, "\n")
obj <- readRDS(input_rds)
meta <- get_metadata_df(obj)

cat("Data loaded. Found", nrow(meta), "cells with", ncol(meta), "metadata columns.\n")

# Detect columns
cluster_col <- choose_first_available(colnames(meta), c("^seurat_clusters$", "cluster$", "cluster_id$"))
if (is.na(cluster_col)) cluster_col <- "seurat_clusters"

perk_col <- choose_first_available(colnames(meta), c("perk", "eif2ak3", "perk_score"))
ire1_col <- choose_first_available(colnames(meta), c("ire1", "ern1", "xbp1", "ire1_score"))
atf6_col <- choose_first_available(colnames(meta), c("atf6", "atf6_score"))

cat("Using columns - Cluster:", cluster_col, "| UPR:", perk_col, ire1_col, atf6_col, "\n")

# Create UPR scores and classify states
meta <- meta %>% mutate(
  UPR_PERK = .data[[perk_col]],
  UPR_IRE1 = .data[[ire1_col]],
  UPR_ATF6 = .data[[atf6_col]],
  cluster_factor = factor(.data[[cluster_col]])
)

meta <- meta %>% mutate(
  UPR_State = classify_rheostat_state(., perk = "UPR_PERK", ire1 = "UPR_IRE1", atf6 = "UPR_ATF6")
)

cat("Creating enhanced Figure 2...\n")

# Calculate correlations using Pearson method
cor_perk_ire1 <- calculate_correlation_with_p(meta$UPR_PERK, meta$UPR_IRE1, method = "pearson")
cor_perk_atf6 <- calculate_correlation_with_p(meta$UPR_PERK, meta$UPR_ATF6, method = "pearson")
cor_ire1_atf6 <- calculate_correlation_with_p(meta$UPR_IRE1, meta$UPR_ATF6, method = "pearson")

# Panel A: PERK vs IRE1 correlation
p2A <- ggplot(meta, aes(UPR_PERK, UPR_IRE1)) +
  stat_density_2d_filled(contour_var = "ndensity", alpha = 0.7) +
  scale_fill_viridis_d(option = "plasma", name = "Density") +
  geom_smooth(method = "lm", se = TRUE, color = "white", fill = "black", 
              alpha = 0.3, linewidth = 1.2) +
  annotate("text", x = 0.05, y = 0.95, 
           label = paste0("r = ", round(cor_perk_ire1$cor, 3), " ", 
                         format_p_value(cor_perk_ire1$p)),
           hjust = 0, vjust = 1, size = 5.5, color = "white", 
           fontface = "bold") +
  labs(title = "PERK - IRE1 Correlation", 
       subtitle = "Inter-arm UPR coordination",
       x = "PERK Activation Score", y = "IRE1 Activation Score") +
  theme_publication_enhanced(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11),
    legend.position = "right"
  )

# Panel B: PERK vs ATF6 correlation
p2B <- ggplot(meta, aes(UPR_PERK, UPR_ATF6)) +
  stat_density_2d_filled(contour_var = "ndensity", alpha = 0.7) +
  scale_fill_viridis_d(option = "plasma", name = "Density") +
  geom_smooth(method = "lm", se = TRUE, color = "white", fill = "grey80", 
              alpha = 0.3, linewidth = 1) +
  annotate("text", x = 0.05, y = max(meta$UPR_ATF6, na.rm = TRUE) * 0.95, 
           label = paste0("r = ", round(cor_perk_atf6$cor, 3), " ", 
                         format_p_value(cor_perk_atf6$p)),
           hjust = 0, vjust = 1, size = 5, color = "white", fontface = "bold") +
  labs(title = "PERK - ATF6 Correlation", 
       subtitle = "Inter-arm UPR coordination",
       x = "PERK Activation Score", y = "ATF6 Activation Score") +
  theme_publication_enhanced(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

# Panel C: IRE1 vs ATF6 correlation
p2C <- ggplot(meta, aes(UPR_IRE1, UPR_ATF6)) +
  stat_density_2d_filled(contour_var = "ndensity", alpha = 0.7) +
  scale_fill_viridis_d(option = "plasma", name = "Density") +
  geom_smooth(method = "lm", se = TRUE, color = "white", fill = "grey80", 
              alpha = 0.3, linewidth = 1) +
  annotate("text", x = 0.05, y = max(meta$UPR_ATF6, na.rm = TRUE) * 0.95, 
           label = paste0("r = ", round(cor_ire1_atf6$cor, 3), " ", 
                         format_p_value(cor_ire1_atf6$p)),
           hjust = 0, vjust = 1, size = 5, color = "white", fontface = "bold") +
  labs(title = "IRE1 - ATF6 Correlation", 
       subtitle = "Inter-arm UPR coordination",
       x = "IRE1 Activation Score", y = "ATF6 Activation Score") +
  theme_publication_enhanced(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    legend.position = "right"
  )

# Panel D: Enrichment heatmap
global_upr_props <- meta %>%
  count(UPR_State) %>%
  mutate(global_prop = n / sum(n)) %>%
  select(UPR_State, global_prop)

upr_heatmap_data <- meta %>%
  count(cluster_factor, UPR_State) %>%
  group_by(cluster_factor) %>%
  mutate(
    total = sum(n),
    fraction = n / total
  ) %>%
  ungroup() %>%
  left_join(global_upr_props, by = "UPR_State") %>%
  mutate(
    enrichment = fraction / global_prop,
    log2_enrichment = log2(enrichment)
  )

p2D <- upr_heatmap_data %>%
  ggplot(aes(x = UPR_State, y = cluster_factor, fill = log2_enrichment)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_distiller(
    palette = "RdBu", 
    direction = 1,
    name = "Enrichment\nScore\n(log2)",
    limits = c(-2, 2),
    oob = scales::squish,
    labels = function(x) round(x, 1)
  ) +
  scale_x_discrete(labels = function(x) gsub("-", "-\n", x)) +
  labs(title = "Cluster-UPR State Enrichment Map", 
       subtitle = "Log2 transformed enrichment scores",
       x = "UPR Rheostat State", y = "Cluster") +
  theme_publication_enhanced(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    panel.grid = element_blank()
  )

# Panel E: Placeholder for LIANA (since no data available) - REMOVED

# Combine panels (without Panel E)
fig2_layout <- (p2A | p2B | p2C) / p2D + 
  plot_annotation(
    tag_levels = list(c('a', 'b', 'c', 'd')),
    title = "UPR Rheostat Hypothesis - Arm Coordination and Cellular Communication",
    subtitle = "Correlations between UPR arms and cluster-specific enrichment patterns",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", color = "grey10", hjust = 0),
      plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0),
      plot.margin = margin(20, 20, 20, 20),
      plot.tag = element_text(size = 18, face = "bold", color = "grey10")
    )
  ) +
  plot_layout(heights = c(1, 1))

# Save figure
ggsave(paste0(output_name, ".pdf"), fig2_layout, 
       width = 16, height = 12, units = "in", dpi = 300)
ggsave(paste0(output_name, ".png"), fig2_layout, 
       width = 16, height = 12, units = "in", dpi = 300, bg = "white")

cat("Enhanced Figure 2 saved\n")

# Export data
correlation_data <- data.frame(
  Comparison = c("PERK-IRE1", "PERK-ATF6", "IRE1-ATF6"),
  Correlation = c(cor_perk_ire1$cor, cor_perk_atf6$cor, cor_ire1_atf6$cor),
  P_value = c(cor_perk_ire1$p, cor_perk_atf6$p, cor_ire1_atf6$p),
  Significance = c(format_p_value(cor_perk_ire1$p), 
                   format_p_value(cor_perk_atf6$p), 
                   format_p_value(cor_ire1_atf6$p)),
  Method = "Pearson"
)

enrichment_data_wide <- upr_heatmap_data %>%
  select(cluster_factor, UPR_State, enrichment, log2_enrichment) %>%
  pivot_wider(
    names_from = UPR_State, 
    values_from = c(enrichment, log2_enrichment),
    names_sep = "_"
  ) %>%
  arrange(as.numeric(as.character(cluster_factor)))

upr_raw_data <- meta %>%
  select(cell_id, cluster_factor, UPR_PERK, UPR_IRE1, UPR_ATF6, UPR_State) %>%
  arrange(as.numeric(as.character(cluster_factor)))

# Export files
write.csv(correlation_data, paste0(output_name, "_correlations.csv"), row.names = FALSE)
write.csv(enrichment_data_wide, paste0(output_name, "_enrichment_heatmap.csv"), row.names = FALSE)
write.csv(upr_heatmap_data, paste0(output_name, "_enrichment_long.csv"), row.names = FALSE)
write.csv(upr_raw_data, paste0(output_name, "_upr_raw_data.csv"), row.names = FALSE)

write.table(correlation_data, paste0(output_name, "_correlations.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(enrichment_data_wide, paste0(output_name, "_enrichment_heatmap.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(upr_raw_data, paste0(output_name, "_upr_raw_data.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Summary statistics
cat("\n=== FIGURE 2 SUMMARY STATISTICS ===\n")
cat("Total cells analyzed:", nrow(meta), "\n")
cat("Number of clusters:", length(unique(meta$cluster_factor)), "\n")

cat("\nUPR Arm Correlations (Pearson):\n")
cat("  PERK-IRE1: r =", round(cor_perk_ire1$cor, 3), format_p_value(cor_perk_ire1$p), "\n")
cat("  PERK-ATF6: r =", round(cor_perk_atf6$cor, 3), format_p_value(cor_perk_atf6$p), "\n")
cat("  IRE1-ATF6: r =", round(cor_ire1_atf6$cor, 3), format_p_value(cor_ire1_atf6$p), "\n")

cat("\nData files exported:\n")
cat("- Correlations: ", output_name, "_correlations.csv/.txt\n", sep = "")
cat("- Enrichment: ", output_name, "_enrichment_heatmap.csv/.txt\n", sep = "")
cat("- Raw data: ", output_name, "_upr_raw_data.csv/.txt\n", sep = "")

cat("\nScript completed successfully!\n")