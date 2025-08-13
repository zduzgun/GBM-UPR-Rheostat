#!/usr/bin/env Rscript
# Complete Figure 1: Glioblastoma Heterogeneity and UPR Rheostat Dynamics
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
  library(ggrepel)
  library(cowplot)
})

# Check for optional libraries
quietly_require <- function(pkg) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) return(FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}

# Input file
input_rds <- "glioblastoma_with_corrected_rheostat.rds"
output_name <- "Figure1_Enhanced_GBM_Heterogeneity"

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

get_embeddings_df <- function(x) {
  if (inherits(x, "Seurat")) {
    red_name <- NULL
    for (nm in c("umap", "UMAP", "tsne", "TSNE", "pca", "PCA")) {
      if (nm %in% names(x@reductions)) { red_name <- nm; break }
    }
    if (is.null(red_name)) {
      # Generate UMAP if missing
      if (quietly_require("Seurat")) {
        DefaultAssay(x) <- Seurat::DefaultAssay(x)
        if (!("pca" %in% names(x@reductions))) x <- Seurat::RunPCA(x, npcs = 30, verbose = FALSE)
        x <- Seurat::RunUMAP(x, reduction = "pca", dims = 1:30, verbose = FALSE)
        red_name <- "umap"
      } else {
        stop("No dimensionality reduction found and Seurat not available to generate UMAP")
      }
    }
    emb <- Seurat::Embeddings(x, reduction = red_name) %>% 
           as.data.frame() %>% 
           tibble::rownames_to_column(var = "cell_id")
    colnames(emb)[2:3] <- c("Dim1", "Dim2")
    attr(emb, "reduction") <- red_name
    return(list(emb = emb, obj = x))
  }
  if (inherits(x, "SingleCellExperiment")) {
    if (!quietly_require("SingleCellExperiment")) stop("SingleCellExperiment required.")
    red_name <- NULL
    for (nm in c("UMAP", "TSNE", "PCA")) {
      if (nm %in% names(SingleCellExperiment::reducedDims(x))) { red_name <- nm; break }
    }
    if (is.null(red_name)) stop("UMAP/TSNE/PCA not found in SCE.")
    emb <- as.data.frame(SingleCellExperiment::reducedDim(x, red_name)) %>% 
           tibble::rownames_to_column(var = "cell_id")
    colnames(emb)[2:3] <- c("Dim1", "Dim2")
    attr(emb, "reduction") <- tolower(red_name)
    return(list(emb = emb, obj = x))
  }
  stop("Object type not supported.")
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

# Enhanced publication theme with modern styling
theme_publication_enhanced <- function(base_size = 12, base_family = "") {
  theme_minimal(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Panel and background
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "grey20", fill = NA, linewidth = 0.8),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92", linewidth = 0.4),
      
      # Text elements
      text = element_text(colour = "grey15"),
      axis.text = element_text(colour = "grey25", size = rel(0.9)),
      axis.title = element_text(colour = "grey15", size = rel(1.0), face = "bold"),
      plot.title = element_text(colour = "grey10", size = rel(1.3), face = "bold", margin = margin(b = 10)),
      plot.subtitle = element_text(colour = "grey40", size = rel(0.9), margin = margin(b = 15)),
      
      # Legend
      legend.key = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.size = unit(0.5, "cm"),
      legend.spacing = unit(0.3, "cm"),
      legend.title = element_text(colour = "grey15", size = rel(0.9), face = "bold"),
      legend.text = element_text(colour = "grey25", size = rel(0.8)),
      legend.margin = margin(l = 10),
      legend.box.margin = margin(0, 0, 0, 0),
      
      # Plot margins and spacing
      plot.margin = unit(c(1.2, 1.2, 1.2, 1.2), "lines"),
      axis.line = element_line(colour = "grey30", linewidth = 0.6),
      axis.ticks = element_line(colour = "grey30", linewidth = 0.5),
      axis.ticks.length = unit(0.15, "cm"),
      
      # Strip for facets
      strip.background = element_rect(colour = "grey20", fill = "grey95", linewidth = 0.6),
      strip.text = element_text(colour = "grey10", size = rel(0.9), face = "bold", margin = margin(4, 4, 4, 4))
    )
}

# Enhanced color palettes
state_palette_enhanced <- c(
  "All-low" = "#E8E8E8",          # Light grey
  "PERK-high" = "#D32F2F",        # Deep red
  "IRE1-high" = "#1976D2",        # Deep blue  
  "ATF6-high" = "#388E3C",        # Deep green
  "Multi-arm-high" = "#7B1FA2",   # Deep purple
  "Intermediate/Mixed" = "#F57C00" # Deep orange
)

# Modern cluster colors with better contrast
cluster_colors_enhanced <- c(
  "#E53E3E", "#3182CE", "#38A169", "#805AD5", "#DD6B20", 
  "#D69E2E", "#319795", "#E53E3E", "#9F7AEA"
)

# Load and process data
cat("Loading data from:", input_rds, "\n")
obj <- readRDS(input_rds)
meta <- get_metadata_df(obj)
emb_list <- get_embeddings_df(obj)
emb <- emb_list$emb
obj <- emb_list$obj

cat("Data loaded. Found", nrow(meta), "cells with", ncol(meta), "metadata columns.\n")
cat("Embeddings type:", attr(emb, "reduction"), "\n")

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

# Merge with coordinates
plot_df <- emb %>% left_join(meta, by = "cell_id")

# Set up cluster palette
n_clusters <- length(unique(plot_df$cluster_factor))
if (n_clusters <= length(cluster_colors_enhanced)) {
  cluster_palette <- setNames(cluster_colors_enhanced[1:n_clusters], sort(unique(plot_df$cluster_factor)))
} else {
  cluster_palette <- setNames(rainbow(n_clusters, s = 0.8, v = 0.8), sort(unique(plot_df$cluster_factor)))
}

cat("Creating enhanced Figure 1...\n")

# ============================================================================
# Panel A: Cluster manifold with enhanced styling
# ============================================================================
p1A <- ggplot(plot_df, aes(Dim1, Dim2, color = cluster_factor)) +
  geom_point(size = 1.0, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = cluster_palette, name = "Cluster") +
  labs(title = "Manifold Map by Clusters", 
       subtitle = paste0("n = ", nrow(plot_df), " cells | ", toupper(attr(emb, "reduction")), " embedding")) +
  theme_void() + 
  theme_publication_enhanced(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

# ============================================================================
# Panel B: UPR states with enhanced visualization
# ============================================================================
p1B <- ggplot(plot_df, aes(Dim1, Dim2, color = UPR_State)) +
  geom_point(size = 1.0, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = state_palette_enhanced, drop = FALSE, name = "UPR State") +
  labs(title = "UPR Rheostat States", 
       subtitle = "Cellular UPR arm activation profiles") +
  theme_void() + 
  theme_publication_enhanced(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1))

# ============================================================================
# Panel C: Enhanced violin plots with better styling
# ============================================================================
long_upr <- plot_df %>%
  select(cell_id, cluster_factor, UPR_PERK, UPR_IRE1, UPR_ATF6) %>%
  pivot_longer(cols = c(UPR_PERK, UPR_IRE1, UPR_ATF6), names_to = "Arm", values_to = "Score") %>%
  mutate(Arm = factor(Arm, levels = c("UPR_PERK", "UPR_IRE1", "UPR_ATF6"),
                      labels = c("PERK", "IRE1", "ATF6")))

# Calculate statistics for annotation
upr_stats <- long_upr %>%
  group_by(cluster_factor, Arm) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p1C <- ggplot(long_upr, aes(x = cluster_factor, y = Score, fill = Arm)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, color = "white", linewidth = 0.3) +
  geom_boxplot(width = 0.12, outlier.size = 0.6, outlier.alpha = 0.6, 
               fill = "white", color = "grey20", linewidth = 0.3, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 21, size = 2.2, 
               fill = "#FFD700", color = "grey20", stroke = 0.6) +
  scale_fill_manual(values = c("PERK" = state_palette_enhanced["PERK-high"], 
                               "IRE1" = state_palette_enhanced["IRE1-high"], 
                               "ATF6" = state_palette_enhanced["ATF6-high"]), 
                    name = "UPR Arm") +
  labs(title = "Cluster-based UPR Arm Activation Profiles", 
       subtitle = "Distribution of UPR pathway activation across cellular clusters",
       x = "Cluster", y = "UPR Arm Activation Score") +
  theme_publication_enhanced(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10),
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    strip.text = element_text(face = "bold", size = 11),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom"
  ) +
  facet_wrap(~Arm, scales = "free_y", nrow = 1) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(alpha = 0.8)))

# ============================================================================
# Panel D: Enhanced stacked bar chart with better styling
# ============================================================================
upr_composition <- plot_df %>%
  count(cluster_factor, UPR_State) %>%
  group_by(cluster_factor) %>% 
  mutate(
    frac = n / sum(n),
    total_cells = sum(n),
    cluster_label = paste0("Cluster ", cluster_factor, "\n(n=", total_cells, ")")
  ) %>% 
  ungroup()

# Calculate cluster order by total cells (largest first)
cluster_order <- upr_composition %>%
  group_by(cluster_factor) %>%
  summarise(total = first(total_cells), .groups = "drop") %>%
  arrange(desc(total)) %>%
  pull(cluster_factor)

upr_composition <- upr_composition %>%
  mutate(cluster_factor = factor(cluster_factor, levels = cluster_order))

p1D <- upr_composition %>%
  ggplot(aes(x = cluster_factor, y = frac, fill = UPR_State)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.4, alpha = 0.9) +
  scale_fill_manual(values = state_palette_enhanced, drop = FALSE, name = "UPR State") +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1), expand = c(0.01, 0.01)) +
  scale_x_discrete(labels = function(x) paste0("Cluster ", x, "\n(n=", 
                                               upr_composition$total_cells[match(x, upr_composition$cluster_factor)], ")")) +
  labs(title = "UPR Rheostat Composition in Cell Clusters", 
       subtitle = "Proportional distribution of UPR activation states",
       x = "Cell Cluster", 
       y = "Proportional Distribution") +
  theme_publication_enhanced(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    axis.text.y = element_text(size = 9),
    axis.text.x = element_text(size = 9),
    axis.title = element_text(size = 10, face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.x = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 0.9)))

# ============================================================================
# Combine panels with enhanced layout
# ============================================================================
fig1_layout <- (p1A | p1B) / p1C / p1D + 
  plot_annotation(
    tag_levels = list(c('a', 'b', 'c', 'd')),
    title = "Glioblastoma Heterogeneity and UPR Rheostat Dynamics",
    subtitle = "Single-cell analysis reveals diverse UPR activation patterns across tumor cell populations",
    theme = theme(
      plot.title = element_text(size = 16, face = "bold", color = "grey10", hjust = 0),
      plot.subtitle = element_text(size = 12, color = "grey40", hjust = 0),
      plot.margin = margin(20, 20, 20, 20),
      plot.tag = element_text(size = 18, face = "bold", color = "grey10")
    ),
    tag_prefix = '',
    tag_suffix = '',
    tag_sep = ''
  ) +
  plot_layout(heights = c(1, 1.3, 1))

# Save enhanced figure
ggsave(paste0(output_name, ".pdf"), fig1_layout, 
       width = 16, height = 18, units = "in", dpi = 300)
ggsave(paste0(output_name, ".png"), fig1_layout, 
       width = 16, height = 18, units = "in", dpi = 300, bg = "white")

# Save high-resolution version
ggsave(paste0(output_name, "_highres.png"), fig1_layout, 
       width = 16, height = 18, units = "in", dpi = 600, bg = "white")

cat("Enhanced Figure 1 saved in multiple formats\n")

# ============================================================================
# Summary statistics and data export
# ============================================================================
cat("\n=== FIGURE 1 SUMMARY STATISTICS ===\n")
cat("Total cells analyzed:", nrow(plot_df), "\n")
cat("Number of clusters:", length(unique(plot_df$cluster_factor)), "\n")
cat("Embedding type:", toupper(attr(emb, "reduction")), "\n")

# Overall UPR distribution
overall_upr <- table(plot_df$UPR_State)
overall_upr_pct <- prop.table(overall_upr) * 100
cat("\nOverall UPR state distribution:\n")
for (i in seq_along(overall_upr_pct)) {
  cat("  ", names(overall_upr_pct)[i], ": ", overall_upr[i], " cells (", 
      round(overall_upr_pct[i], 1), "%)\n", sep = "")
}

# Cluster composition summary
cat("\nCluster composition:\n")
cluster_summary <- plot_df %>%
  count(cluster_factor, name = "cells") %>%
  mutate(percentage = round(cells / sum(cells) * 100, 1)) %>%
  arrange(desc(cells))
print(cluster_summary)

# ============================================================================
# Export comprehensive data tables (matching your original calc_upr_percent.R output)
# ============================================================================

# Wide format composition table (matching upr_state_percentages.txt)
composition_table_wide <- upr_composition %>%
  select(cluster_factor, UPR_State, frac) %>%
  mutate(percentage = round(frac * 100, 2)) %>%
  pivot_wider(names_from = UPR_State, values_from = percentage, values_fill = 0) %>%
  left_join(cluster_summary %>% select(cluster_factor, Total_Cells = cells), by = "cluster_factor") %>%
  relocate(Total_Cells, .after = cluster_factor) %>%
  arrange(as.numeric(as.character(cluster_factor)))

# Long format composition table (matching upr_state_percentages_long_format.txt)
composition_table_long <- upr_composition %>%
  mutate(percentage = round(frac * 100, 15)) %>%  # High precision to match original
  select(seurat_clusters = cluster_factor, UPR_State, n, total_cells, frac, percentage) %>%
  arrange(as.numeric(as.character(seurat_clusters)), UPR_State)

# Cluster-based table (matching upr_state_percentages_by_cluster.txt)
composition_table_cluster <- composition_table_wide %>%
  rename(seurat_clusters = cluster_factor)

# Print the tables in console (matching your original output)
cat("\n=== UPR STATE PERCENTAGES BY CLUSTER (Wide Format) ===\n")
print(composition_table_wide)

cat("\n=== UPR STATE PERCENTAGES (Long Format) - First 20 rows ===\n")
print(head(composition_table_long, 20))

cat("\n=== UPR STATE PERCENTAGES BY CLUSTER ===\n")
print(composition_table_cluster)

# Export all data tables
write.table(composition_table_wide, paste0(output_name, "_composition_wide.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(composition_table_long, paste0(output_name, "_composition_long.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(composition_table_cluster, paste0(output_name, "_composition_by_cluster.txt"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# Also export as CSV for convenience
write.csv(composition_table_wide, paste0(output_name, "_composition_wide.csv"), row.names = FALSE)
write.csv(composition_table_long, paste0(output_name, "_composition_long.csv"), row.names = FALSE)
write.csv(upr_stats, paste0(output_name, "_violin_statistics.csv"), row.names = FALSE)

cat("\nData tables exported:\n")
cat("TXT format (tab-separated, matching original calc_upr_percent.R):\n")
cat("- ", output_name, "_composition_wide.txt\n", sep = "")
cat("- ", output_name, "_composition_long.txt\n", sep = "")
cat("- ", output_name, "_composition_by_cluster.txt\n", sep = "")
cat("\nCSV format (for Excel/analysis):\n")
cat("- ", output_name, "_composition_wide.csv\n", sep = "")
cat("- ", output_name, "_composition_long.csv\n", sep = "")
cat("- ", output_name, "_violin_statistics.csv\n", sep = "")

cat("\nScript completed successfully!\n")
cat("Files generated:\n")
cat("- ", output_name, ".pdf (standard resolution)\n", sep = "")
cat("- ", output_name, ".png (standard resolution)\n", sep = "")
cat("- ", output_name, "_highres.png (high resolution)\n", sep = "")