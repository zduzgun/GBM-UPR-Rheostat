#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
  library(ggridges)
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(ggrepel)
  library(corrplot)
})

quietly_require <- function(pkg) {
  if (!suppressWarnings(requireNamespace(pkg, quietly = TRUE))) return(FALSE)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  TRUE
}

log_message <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste0(..., collapse = " ")))
  flush.console()
}

stop_with_hint <- function(msg) {
  log_message(msg)
  quit(status = 1)
}

# Statistical test functions
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

# Enhanced publication theme function with better contrast
theme_publication <- function(base_size = 12, base_family = "") {
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      panel.background = element_rect(fill = "white", colour = NA),
      plot.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
      strip.background = element_rect(colour = "black", fill = "grey90", linewidth = 0.8),
      strip.text = element_text(colour = "black", size = rel(0.9), face = "bold"),
      axis.text = element_text(colour = "black", size = rel(0.9)),
      axis.title = element_text(colour = "black", size = rel(1.0), face = "bold"),
      plot.title = element_text(colour = "black", size = rel(1.2), face = "bold"),
      plot.subtitle = element_text(colour = "grey20", size = rel(0.9)),
      legend.key = element_rect(fill = "white", colour = NA),
      legend.position = "right",
      legend.direction = "vertical",
      legend.key.size = unit(0.4, "cm"),
      legend.spacing = unit(0.2, "cm"),
      legend.title = element_text(colour = "black", size = rel(0.9), face = "bold"),
      legend.text = element_text(colour = "black", size = rel(0.8)),
      legend.margin = margin(5, 5, 5, 5),
      legend.box.margin = margin(0, 0, 0, 0),
      plot.margin = unit(c(1, 1, 1, 1), "lines"),
      axis.line = element_line(colour = "black", linewidth = 0.8),
      axis.ticks = element_line(colour = "black", linewidth = 0.6),
      axis.ticks.length = unit(0.2, "cm")
    )
}

# ----------------------
# Input files and output names
# ----------------------
input_rds <- "glioblastoma_with_corrected_reostat.rds"

out_fig1 <- "Figure1_GBM_Heterogeneity.pdf"
out_fig2 <- "Figure2_UPR_Reostat_Hypothesis.pdf"
out_fig3 <- "Figure3_UPR_Stress_Integration_combined.pdf"

# LIANA-derived possible files (will be used in Figure 2 if available)
liana_top_file <- "liana_top100_interactions.csv"
cluster_sender_file <- "cluster_sender_summary.csv"
cluster_receiver_file <- "cluster_receiver_summary.csv"

if (!file.exists(input_rds)) {
  stop_with_hint(sprintf("Expected RDS file not found: %s", input_rds))
}

log_message("Loading object:", input_rds)
obj <- readRDS(input_rds)

# ----------------------
# Object type-dependent helpers
# ----------------------
is_seurat <- inherits(obj, "Seurat")
is_sce <- inherits(obj, "SingleCellExperiment")

get_metadata_df <- function(x) {
  if (inherits(x, "Seurat")) return(x@meta.data %>% tibble::rownames_to_column(var = "cell_id"))
  if (inherits(x, "SingleCellExperiment")) return(as.data.frame(SummarizedExperiment::colData(x)) %>% tibble::rownames_to_column(var = "cell_id"))
  stop("Object type not supported. Expected Seurat or SingleCellExperiment.")
}

get_embeddings_df <- function(x) {
  if (inherits(x, "Seurat")) {
    red_name <- NULL
    for (nm in c("umap", "UMAP", "tsne", "pca")) {
      if (nm %in% names(x@reductions)) { red_name <- nm; break }
    }
    if (is.null(red_name)) {
      # If UMAP is missing, try to quickly generate PCA->UMAP
      if (!quietly_require("Seurat")) stop("Seurat required.")
      DefaultAssay(x) <- Seurat::DefaultAssay(x)
      if (!("pca" %in% names(x@reductions))) x <- Seurat::RunPCA(x, npcs = 30, verbose = FALSE)
      x <- Seurat::RunUMAP(x, reduction = "pca", dims = 1:30, verbose = FALSE)
      red_name <- "umap"
    }
    emb <- Seurat::Embeddings(x, reduction = red_name) %>% as.data.frame() %>% tibble::rownames_to_column(var = "cell_id")
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
    if (is.null(red_name)) stop("UMAP/TSNE/PCA not found in SCE. Please add reduced dimensions.")
    emb <- as.data.frame(SingleCellExperiment::reducedDim(x, red_name)) %>% tibble::rownames_to_column(var = "cell_id")
    colnames(emb)[2:3] <- c("Dim1", "Dim2")
    attr(emb, "reduction") <- tolower(red_name)
    return(list(emb = emb, obj = x))
  }
  stop("Object type not supported.")
}

choose_first_available <- function(cn, patterns) {
  for (pt in patterns) {
    idx <- which(grepl(pt, cn, ignore.case = TRUE))
    if (length(idx) > 0) return(cn[idx[1]])
  }
  return(NA_character_)
}

# ----------------------
# Metadata and coordinates
# ----------------------
meta <- get_metadata_df(obj)
emb_list <- get_embeddings_df(obj)
emb <- emb_list$emb
obj <- emb_list$obj

# Try to capture cluster and cell type columns
cluster_col <- choose_first_available(colnames(meta), c("^seurat_clusters$", "cluster$", "cluster_id$", "snn_res", "leiden", "louvain"))
if (is.na(cluster_col)) cluster_col <- choose_first_available(colnames(meta), c("cluster", "group", "subcluster"))
if (is.na(cluster_col)) cluster_col <- "cluster"

celltype_col <- choose_first_available(colnames(meta), c("celltype", "cell_type", "CellType", "annotation", "cell_identity", "cell_class"))
if (is.na(celltype_col)) celltype_col <- cluster_col

# Capture UPR arms and total/reostat approaches
perk_col <- choose_first_available(colnames(meta), c("perk", "eif2ak3", "perk_score"))
ire1_col <- choose_first_available(colnames(meta), c("ire1", "ern1", "xbp1", "ire1_score"))
atf6_col <- choose_first_available(colnames(meta), c("atf6", "atf6_score"))
reostat_col <- choose_first_available(colnames(meta), c("reostat", "upr_reostat", "upr_reostat_score"))
upr_total_col <- choose_first_available(colnames(meta), c("upr_total", "upr_score", "unfolded_protein_response"))

if (is.na(perk_col) || is.na(ire1_col) || is.na(atf6_col)) {
  stop_with_hint("UPR arms (PERK, IRE1, ATF6) not found in metadata. Make sure relevant score columns are present in 'glioblastoma_with_corrected_reostat.rds'.")
}

# Reostat state classification
classify_reostat_state <- function(df, perk, ire1, atf6, q = 0.75, low_q = 0.25) {
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

meta <- meta %>% mutate(
  UPR_PERK = .data[[perk_col]],
  UPR_IRE1 = .data[[ire1_col]],
  UPR_ATF6 = .data[[atf6_col]]
)

meta <- meta %>% mutate(
  UPR_State = classify_reostat_state(., perk = "UPR_PERK", ire1 = "UPR_IRE1", atf6 = "UPR_ATF6")
)

if (!is.na(upr_total_col)) meta$UPR_Total <- meta[[upr_total_col]]
if (!is.na(reostat_col)) meta$UPR_Reostat <- meta[[reostat_col]]

# Merge with UMAP/TSNE coordinates
plot_df <- emb %>% left_join(meta, by = "cell_id")

# Enhanced color palettes (colorblind-friendly) with better contrast
state_palette <- c(
  "All-low" = "#CCCCCC",
  "PERK-high" = "#E31A1C", 
  "IRE1-high" = "#1F78B4",
  "ATF6-high" = "#33A02C",
  "Multi-arm-high" = "#9A3192",
  "Intermediate/Mixed" = "#FF7F00"
)

# Nature-style color palette with higher saturation
nature_colors <- c("#E31A1C", "#1F78B4", "#33A02C", "#6A3D9A", "#FF7F00", 
                   "#B15928", "#A6CEE3", "#FB9A99", "#CAB2D6", "#FDBF6F")

# Additional high-contrast palette for special cases
high_contrast_colors <- c("#000000", "#E31A1C", "#1F78B4", "#33A02C", 
                         "#6A3D9A", "#FF7F00", "#B15928", "#FFFF99")

cluster_palette <- NULL
if (!is.na(cluster_col)) {
  cls <- sort(unique(plot_df[[cluster_col]]))
  n_cls <- length(cls)
  if (n_cls <= length(nature_colors)) {
    cluster_palette <- setNames(nature_colors[1:n_cls], cls)
  } else {
    cluster_palette <- setNames(rainbow(n_cls, s = 0.8, v = 0.8), cls)
  }
}

celltype_palette <- NULL
if (!is.na(celltype_col)) {
  cts <- sort(unique(plot_df[[celltype_col]]))
  n_cts <- length(cts)
  if (n_cts <= length(nature_colors)) {
    celltype_palette <- setNames(nature_colors[1:n_cts], cts)
  } else {
    celltype_palette <- setNames(rainbow(n_cts, s = 0.8, v = 0.8), cts)
  }
}

# ----------------------
# Figure 1: Heterogeneity and UPR states
# ----------------------
log_message("Creating Figure 1...")

# Panel A: Cell type manifold (enhanced visibility)
p1A <- ggplot(plot_df, aes(Dim1, Dim2, color = .data[[celltype_col]])) +
  geom_point(size = 0.8, alpha = 0.8, stroke = 0.1) +
  scale_color_manual(values = celltype_palette) +
  labs(title = "Manifold Map by Cell Types", 
       subtitle = paste0("n = ", nrow(plot_df), " cells"),
       color = "Cell Type") +
  theme_void() + 
  theme_publication(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, stroke = 0)))

# Panel B: UPR reostat states (enhanced visualization)
p1B <- ggplot(plot_df, aes(Dim1, Dim2, color = UPR_State)) +
  geom_point(size = 0.8, alpha = 0.9, stroke = 0.1) +
  scale_color_manual(values = state_palette, drop = FALSE) +
  labs(title = "UPR Reostat States", 
       subtitle = "Cellular UPR arm activation profiles",
       color = "UPR State") +
  theme_void() + 
  theme_publication(base_size = 12) +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14, color = "black"),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.key.size = unit(0.6, "cm"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, stroke = 0)))

# Panel C: UPR arm distributions (enhanced violin plot)
long_upr <- plot_df %>%
  select(cell_id, !!sym(cluster_col), UPR_PERK, UPR_IRE1, UPR_ATF6) %>%
  pivot_longer(cols = c(UPR_PERK, UPR_IRE1, UPR_ATF6), names_to = "Arm", values_to = "Score") %>%
  mutate(Arm = factor(Arm, levels = c("UPR_PERK", "UPR_IRE1", "UPR_ATF6"),
                      labels = c("PERK", "IRE1", "ATF6")))

# For statistical comparisons
upr_stats <- long_upr %>%
  group_by(!!sym(cluster_col), Arm) %>%
  summarise(
    mean_score = mean(Score, na.rm = TRUE),
    median_score = median(Score, na.rm = TRUE),
    .groups = "drop"
  )

p1C <- ggplot(long_upr, aes(x = .data[[cluster_col]], y = Score, fill = Arm)) +
  geom_violin(scale = "width", trim = TRUE, alpha = 0.7, color = "white", linewidth = 0.4) +
  geom_boxplot(width = 0.15, outlier.size = 0.8, outlier.alpha = 0.7, 
               fill = "white", color = "black", linewidth = 0.4) +
  stat_summary(fun = mean, geom = "point", shape = 23, size = 2.5, 
               fill = "yellow", color = "black", stroke = 0.8) +
  scale_fill_manual(values = c("PERK" = state_palette["PERK-high"], 
                               "IRE1" = state_palette["IRE1-high"], 
                               "ATF6" = state_palette["ATF6-high"])) +
  labs(title = "Cluster-based UPR Arm Activation Profiles", 
       subtitle = "Colorful violin plots show distribution density per cluster",
       x = "Cluster", y = "UPR Arm Activation Score", fill = "UPR Arm") +
  theme_publication(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    strip.text = element_text(face = "bold", size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom"
  ) +
  facet_wrap(~Arm, scales = "free_y", nrow = 1)

# Panel D: UPR state distribution per cell type (enhanced stacked bar)
upr_composition <- plot_df %>%
  count(!!sym(celltype_col), UPR_State) %>%
  group_by(!!sym(celltype_col)) %>% 
  mutate(
    frac = n / sum(n),
    total_cells = sum(n),
    celltype_label = paste0(.data[[celltype_col]], "\n(n=", total_cells, ")")
  ) %>% 
  ungroup()

p1D <- upr_composition %>%
  ggplot(aes(x = reorder(celltype_label, total_cells), y = frac, fill = UPR_State)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.3, alpha = 0.9) +
  scale_fill_manual(values = state_palette, drop = FALSE) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(), expand = c(0, 0)) +
  labs(title = "UPR Reostat Composition in Cell Types", 
       subtitle = "Proportional distribution of UPR states in each cell type",
       x = "Cell Type (cell count)", y = "Proportional Distribution", fill = "UPR State") +
  theme_publication(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey50"),
    axis.text.y = element_text(size = 8)
  ) +
  guides(fill = guide_legend(ncol = 2))

# Enhanced figure layout
fig1 <- (p1A | p1B) / (p1C) / (p1D) + 
  plot_annotation(
    tag_levels = 'A',
    title = "Figure 1: Glioblastoma Heterogeneity and UPR Reostat Dynamics",
    subtitle = "Cellular and spatial distribution of UPR arm activation in single-cell RNA sequencing data",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "grey30")
    )
  ) +
  plot_layout(heights = c(1, 1.2, 0.8))

# Save as PDF
ggsave(out_fig1, fig1, width = 14, height = 16, units = "in", dpi = 300)
log_message("Figure 1 PDF saved:", out_fig1)

# Save as high-resolution PNG
png_fig1 <- gsub("\\.pdf$", "_1200dpi.png", out_fig1)
ggsave(png_fig1, fig1, width = 14, height = 16, units = "in", dpi = 1200, bg = "white")
log_message("Figure 1 PNG saved:", png_fig1)

# ----------------------
# Figure 2: UPR reostat hypothesis - arm relationships and communication context
# ----------------------
log_message("Creating Figure 2...")

# Calculate statistical correlations
cor_perk_ire1 <- calculate_correlation_with_p(plot_df$UPR_PERK, plot_df$UPR_IRE1)
cor_perk_atf6 <- calculate_correlation_with_p(plot_df$UPR_PERK, plot_df$UPR_ATF6)
cor_ire1_atf6 <- calculate_correlation_with_p(plot_df$UPR_IRE1, plot_df$UPR_ATF6)

use_hexbin <- requireNamespace("hexbin", quietly = TRUE)

dens_layer <- function() {
  if (use_hexbin) return(ggplot2::stat_binhex(bins = 50))
  ggplot2::stat_density_2d_filled(contour_var = "ndensity", alpha = 0.8)
}

# Panel A: PERK vs IRE1 correlation (enhanced visualization)
p2A <- ggplot(plot_df, aes(UPR_PERK, UPR_IRE1)) +
  dens_layer() + 
  scale_fill_viridis_c(option = "plasma", name = "Density", 
                       guide = guide_colorbar(barwidth = 1, barheight = 8)) +
  geom_smooth(method = "lm", se = TRUE, color = "white", fill = "black", 
              alpha = 0.2, linewidth = 1.5) +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("r = ", round(cor_perk_ire1$cor, 3), " ", 
                         format_p_value(cor_perk_ire1$p)),
           hjust = 1.05, vjust = 1.5, size = 5, color = "white", 
           fontface = "bold") +
  labs(title = "PERK - IRE1 Correlation", 
       subtitle = "Inter-arm UPR coordination",
       x = "PERK Activation Score", y = "IRE1 Activation Score") +
  theme_publication(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 11)
  )

# Panel B: PERK vs ATF6 correlation
p2B <- ggplot(plot_df, aes(UPR_PERK, UPR_ATF6)) +
  dens_layer() + 
  scale_fill_viridis_c(option = "plasma", name = "Density") +
  geom_smooth(method = "lm", se = TRUE, color = "white", fill = "grey80", 
              alpha = 0.3, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("r = ", round(cor_perk_atf6$cor, 3), " ", 
                         format_p_value(cor_perk_atf6$p)),
           hjust = 1.05, vjust = 1.5, size = 5, color = "white", fontface = "bold") +
  labs(title = "PERK - ATF6 Correlation", 
       subtitle = "Inter-arm UPR coordination",
       x = "PERK Activation Score", y = "ATF6 Activation Score") +
  theme_publication(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey50")
  )

# Panel C: IRE1 vs ATF6 correlation
p2C <- ggplot(plot_df, aes(UPR_IRE1, UPR_ATF6)) +
  dens_layer() + 
  scale_fill_viridis_c(option = "plasma", name = "Density") +
  geom_smooth(method = "lm", se = TRUE, color = "white", fill = "grey80", 
              alpha = 0.3, linewidth = 1) +
  annotate("text", x = Inf, y = Inf, 
           label = paste0("r = ", round(cor_ire1_atf6$cor, 3), " ", 
                         format_p_value(cor_ire1_atf6$p)),
           hjust = 1.05, vjust = 1.5, size = 5, color = "white", fontface = "bold") +
  labs(title = "IRE1 - ATF6 Correlation", 
       subtitle = "Inter-arm UPR coordination",
       x = "IRE1 Activation Score", y = "ATF6 Activation Score") +
  theme_publication(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey50")
  )

# Panel D: Heatmap visualization of UPR reostat states
# Calculate global UPR state proportions first
global_upr_props <- plot_df %>%
  count(UPR_State) %>%
  mutate(global_prop = n / sum(n)) %>%
  select(UPR_State, global_prop)

upr_heatmap_data <- plot_df %>%
  group_by(!!sym(cluster_col), UPR_State) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(!!sym(cluster_col)) %>%
  mutate(
    total = sum(count),
    fraction = count / total
  ) %>%
  ungroup() %>%
  left_join(global_upr_props, by = "UPR_State") %>%
  mutate(enrichment = fraction / global_prop)

p2D <- upr_heatmap_data %>%
  ggplot(aes(x = UPR_State, y = .data[[cluster_col]], fill = enrichment)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_distiller(
    palette = "RdBu", 
    direction = 1,
    name = "Enrichment\nScore",
    trans = "log2",
    labels = function(x) round(x, 2)
  ) +
  labs(title = "Cluster-UPR State Enrichment Map", 
       subtitle = "Log2 transformed enrichment scores",
       x = "UPR Reostat State", y = "Cluster") +
  theme_publication(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 11),
    plot.subtitle = element_text(size = 9, color = "grey50")
  )

# Panel E: LIANA ligand-receptor interactions (enhanced visualization)
p2E <- NULL
if (file.exists(liana_top_file)) {
  li_top <- suppressMessages(readr::read_csv(liana_top_file, show_col_types = FALSE))
  score_col <- choose_first_available(colnames(li_top), c("score", "magnitude", "lrscore", "prob", "pvalue", "rank"))
  lig_col <- choose_first_available(colnames(li_top), c("ligand", "ligand_complex", "ligand_symbol"))
  rec_col <- choose_first_available(colnames(li_top), c("receptor", "receptor_complex", "receptor_symbol"))
  source_col <- choose_first_available(colnames(li_top), c("source", "sender", "ligand_cluster"))
  target_col <- choose_first_available(colnames(li_top), c("target", "receiver", "receptor_cluster"))
  
  if (!is.na(score_col) && !is.na(lig_col) && !is.na(rec_col)) {
    li_plot_df <- li_top %>% 
      mutate(
        lr = paste(.data[[lig_col]], .data[[rec_col]], sep = " → "),
        interaction_type = case_when(
          grepl("stress|upr|ER|unfold", paste(.data[[lig_col]], .data[[rec_col]]), ignore.case = TRUE) ~ "Stress-related",
          grepl("growth|prolif|cycle", paste(.data[[lig_col]], .data[[rec_col]]), ignore.case = TRUE) ~ "Growth-related", 
          TRUE ~ "Other"
        )
      ) %>%
      arrange(desc(.data[[score_col]])) %>% 
      slice_head(n = 15) %>%
      mutate(lr = factor(lr, levels = rev(unique(lr))))
    
    p2E <- ggplot(li_plot_df, aes(x = lr, y = .data[[score_col]], fill = interaction_type)) +
      geom_col(alpha = 0.8, color = "white", linewidth = 0.3) +
      scale_fill_manual(values = c("Stress-related" = "#D55E00", 
                                   "Growth-related" = "#0173B2", 
                                   "Other" = "#666666")) +
      coord_flip() +
      labs(title = "Top Ligand-Receptor Interactions", 
           subtitle = "Intercellular communication network",
           x = "Ligand → Receptor", y = paste("Interaction Score (", score_col, ")"),
           fill = "Interaction Type") +
      theme_publication(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 9, color = "grey50"),
        axis.text.y = element_text(size = 8)
      )
  }
}

if (is.null(p2E)) {
  p2E <- ggplot() + 
    theme_void() + 
    labs(title = "LIANA Results Not Found",
         subtitle = "Ligand-receptor analysis not available") +
    theme_publication(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
      plot.subtitle = element_text(size = 9, color = "grey50", hjust = 0.5)
    )
}

# Enhanced Figure 2 layout
fig2 <- (p2A | p2B | p2C) / (p2D | p2E) + 
  plot_annotation(
    tag_levels = 'A',
    title = "Figure 2: UPR Reostat Hypothesis - Arm Coordination and Cellular Communication",
    subtitle = "Correlations between UPR arms and intercellular interaction networks",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "grey30")
    )
  ) +
  plot_layout(heights = c(1, 1))

# Save as PDF
ggsave(out_fig2, fig2, width = 16, height = 12, units = "in", dpi = 300)
log_message("Figure 2 PDF saved:", out_fig2)

# Save as high-resolution PNG
png_fig2 <- gsub("\\.pdf$", "_1200dpi.png", out_fig2)
ggsave(png_fig2, fig2, width = 16, height = 12, units = "in", dpi = 1200, bg = "white")
log_message("Figure 2 PNG saved:", png_fig2)

# ----------------------
# Figure 3: UPR-stress integration and trajectory context
# ----------------------
log_message("Creating Figure 3...")

# Enhanced stress module detection
stress_like_cols <- grep("hypoxia|oxidative|stress|emt|apoptosis|dna_damage|unfolded|proteasome|inflammation|nfkb|tnfa|ifn|il6|glycolysis|oxidative_phosphorylation|cholesterol|lipid|autophagy|angiogenesis|stemness|differentiation",
                         colnames(meta), ignore.case = TRUE, value = TRUE)

sel_cols <- unique(c("UPR_PERK", "UPR_IRE1", "UPR_ATF6", intersect(stress_like_cols, colnames(meta))))

# Panel A: Enhanced correlation matrix
cor_mat_plot <- function(df, cols) {
  if (length(cols) < 3) {
    return(ggplot() + 
           theme_void() + 
           labs(title = "Insufficient Stress Module Metadata",
                subtitle = "Additional stress pathway scores required") +
           theme_publication(base_size = 12) +
           theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
                 plot.subtitle = element_text(color = "grey30", hjust = 0.5, size = 12)))
  }
  
  # Ensure we have valid numeric data
  sub <- df %>% 
    select(all_of(cols)) %>% 
    mutate(across(everything(), as.numeric)) %>%
    select(where(~ !all(is.na(.x)))) # Remove columns with all NAs
  
  if (ncol(sub) < 2) {
    return(ggplot() + theme_void() + labs(title = "Insufficient valid data for correlation"))
  }
  
  cm <- suppressWarnings(cor(sub, use = "pairwise.complete.obs", method = "pearson"))
  
  # Calculate p-values safely
  n_vars <- ncol(sub)
  p_matrix <- matrix(1, n_vars, n_vars) # Initialize with 1s
  for(i in 1:n_vars) {
    for(j in 1:n_vars) {
      if(i != j) {
        tryCatch({
          cor_test <- suppressWarnings(cor.test(sub[[i]], sub[[j]], method = "pearson"))
          if (!is.na(cor_test$p.value)) {
            p_matrix[i,j] <- cor_test$p.value
          }
        }, error = function(e) {
          p_matrix[i,j] <<- 1
        })
      }
    }
  }
  
  # Convert correlation matrix to dataframe
  cm_df <- as.data.frame(as.table(cm))
  colnames(cm_df) <- c("Var1", "Var2", "correlation")
  
  # Add p-values
  p_df <- as.data.frame(as.table(p_matrix))
  colnames(p_df) <- c("Var1", "Var2", "p_value")
  
  plot_df <- cm_df %>% 
    left_join(p_df, by = c("Var1", "Var2")) %>%
    mutate(
      significance = format_p_value(p_value),
      correlation_text = ifelse(Var1 == Var2, "", 
                               paste0(round(correlation, 2), "\n", significance)),
      text_color = ifelse(abs(correlation) > 0.5, "white", "black")
    )
  
  ggplot(plot_df, aes(Var1, Var2, fill = correlation)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = correlation_text, color = I(text_color)), 
              size = 3.5, fontface = "bold") +
    scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), oob = squish, direction = 1,
                         name = "Pearson\nr",
                         guide = guide_colorbar(barwidth = 1.5, barheight = 10)) +
    coord_equal() +
    labs(title = "UPR and Stress Pathway Correlation Matrix", 
         subtitle = "Pearson correlations and significance levels",
         x = NULL, y = NULL) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30")
    )
}

p3A <- cor_mat_plot(meta, sel_cols)

# Panel B: Pseudotime analysis (enhanced trajectory visualization)
pseudotime_col <- choose_first_available(colnames(meta), c("pseudotime", "monocle_pseudotime", "trajectory", "time", "latent_time"))
p3B <- NULL
if (!is.na(pseudotime_col)) {
  # Divide pseudotime into bins
  trajectory_data <- meta %>%
    mutate(
      pseudotime_bin = cut(.data[[pseudotime_col]], breaks = 20, labels = FALSE),
      pseudotime_value = .data[[pseudotime_col]]
    ) %>%
    select(pseudotime_bin, pseudotime_value, UPR_PERK, UPR_IRE1, UPR_ATF6) %>%
    pivot_longer(cols = c(UPR_PERK, UPR_IRE1, UPR_ATF6), names_to = "Arm", values_to = "Score") %>%
    mutate(Arm = factor(Arm, levels = c("UPR_PERK", "UPR_IRE1", "UPR_ATF6"),
                        labels = c("PERK", "IRE1", "ATF6")))
  
  # Calculate mean and confidence interval for each pseudotime bin
  trajectory_summary <- trajectory_data %>%
    group_by(pseudotime_bin, Arm) %>%
    summarise(
      mean_pseudotime = mean(pseudotime_value, na.rm = TRUE),
      mean_score = mean(Score, na.rm = TRUE),
      se_score = sd(Score, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    )
  
  p3B <- ggplot(trajectory_data, aes(x = pseudotime_value, y = Score, color = Arm, fill = Arm)) +
    geom_point(alpha = 0.1, size = 0.2) +
    geom_smooth(se = TRUE, method = "loess", span = 0.3, alpha = 0.2, linewidth = 1.2) +
    scale_color_manual(values = c("PERK" = state_palette["PERK-high"], 
                                  "IRE1" = state_palette["IRE1-high"], 
                                  "ATF6" = state_palette["ATF6-high"])) +
    scale_fill_manual(values = c("PERK" = state_palette["PERK-high"], 
                                 "IRE1" = state_palette["IRE1-high"], 
                                 "ATF6" = state_palette["ATF6-high"])) +
    labs(title = "UPR Dynamics Along Differentiation Trajectory", 
         subtitle = "Trend analysis with LOESS smoothing",
         x = paste("Pseudotime (", pseudotime_col, ")"), 
         y = "UPR Arm Activation Score", 
         color = "UPR Arm", fill = "UPR Arm") +
    theme_publication(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11),
      plot.subtitle = element_text(size = 9, color = "grey50")
    ) +
    facet_wrap(~Arm, scales = "free_y", nrow = 1)
} else {
  # Alternative: UPR arm balance plot
  upr_balance <- plot_df %>%
    mutate(
      PERK_IRE1_ratio = (UPR_PERK + 0.001) / (UPR_IRE1 + 0.001),
      PERK_ATF6_ratio = (UPR_PERK + 0.001) / (UPR_ATF6 + 0.001),
      IRE1_ATF6_ratio = (UPR_IRE1 + 0.001) / (UPR_ATF6 + 0.001)
    ) %>%
    select(!!sym(cluster_col), PERK_IRE1_ratio, PERK_ATF6_ratio, IRE1_ATF6_ratio) %>%
    pivot_longer(cols = c(PERK_IRE1_ratio, PERK_ATF6_ratio, IRE1_ATF6_ratio), 
                 names_to = "Ratio", values_to = "Value") %>%
    mutate(Ratio = case_when(
      Ratio == "PERK_IRE1_ratio" ~ "PERK/IRE1",
      Ratio == "PERK_ATF6_ratio" ~ "PERK/ATF6", 
      Ratio == "IRE1_ATF6_ratio" ~ "IRE1/ATF6"
    )) %>%
    filter(is.finite(Value) & Value > 0)  # Remove infinite, NaN, and non-positive values
  
  p3B <- ggplot(upr_balance, aes(x = .data[[cluster_col]], y = Value, fill = Ratio)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
    scale_fill_manual(values = c("PERK/IRE1" = "#E31A1C", 
                                 "PERK/ATF6" = "#1F78B4", 
                                 "IRE1/ATF6" = "#33A02C"),
                      guide = guide_legend(title = "Arm Ratio")) +
    scale_y_log10(labels = scales::number_format(accuracy = 0.01)) +
    labs(title = "UPR Arm Balance Ratios by Cluster", 
         subtitle = "Log-scale ratios show relative arm dominance",
         x = "Cluster", y = "UPR Arm Ratio (log scale)", fill = "Arm Ratio") +
    theme_publication(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30")
    ) +
    facet_wrap(~Ratio, scales = "free_y", nrow = 1)
}

# Panel C: UPR reostat score manifold visualization
color_col <- if (!is.null(meta$UPR_Total)) "UPR_Total" else if (!is.null(meta$UPR_Reostat)) "UPR_Reostat" else "UPR_State"

if (color_col == "UPR_State") {
  p3C <- ggplot(plot_df, aes(Dim1, Dim2, color = UPR_State)) +
    geom_point(size = 1.2, alpha = 0.9, stroke = 0.1) +
    scale_color_manual(values = state_palette, drop = FALSE) +
    labs(title = "Spatial Distribution of UPR Reostat States", 
         subtitle = "Manifold map of cellular UPR activation profiles",
         color = "UPR State") +
    theme_void() + 
    theme_publication(base_size = 12) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30"),
      legend.key.size = unit(0.6, "cm"),
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1, stroke = 0)))
} else {
  p3C <- ggplot(plot_df, aes(Dim1, Dim2, color = .data[[color_col]])) +
    geom_point(size = 1.2, alpha = 0.9, stroke = 0.1) +
    scale_color_viridis_c(option = "inferno", name = gsub("_", " ", color_col),
                          guide = guide_colorbar(barwidth = 1.5, barheight = 10)) +
    labs(title = paste0("Manifold View of ", gsub("_", " ", color_col)), 
         subtitle = "Continuous UPR activation scores",
         color = gsub("_", " ", color_col)) +
    theme_void() + 
    theme_publication(base_size = 12) +
    theme(
      axis.line = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11, color = "grey30"),
      legend.title = element_text(size = 12, face = "bold")
    )
}

# Panel D: Numerical summary of UPR states
upr_summary_stats <- plot_df %>%
  group_by(UPR_State) %>%
  summarise(
    count = n(),
    mean_perk = mean(UPR_PERK, na.rm = TRUE),
    mean_ire1 = mean(UPR_IRE1, na.rm = TRUE),
    mean_atf6 = mean(UPR_ATF6, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    percentage = count / sum(count) * 100,
    state_label = paste0(UPR_State, "\n(", round(percentage, 1), "%)")
  ) %>%
  pivot_longer(cols = c(mean_perk, mean_ire1, mean_atf6), 
               names_to = "arm", values_to = "mean_score") %>%
  mutate(arm = case_when(
    arm == "mean_perk" ~ "PERK",
    arm == "mean_ire1" ~ "IRE1", 
    arm == "mean_atf6" ~ "ATF6"
  ))

p3D <- ggplot(upr_summary_stats, aes(x = UPR_State, y = mean_score, fill = arm)) +
  geom_col(position = "dodge", alpha = 0.9, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("PERK" = state_palette["PERK-high"], 
                               "IRE1" = state_palette["IRE1-high"], 
                               "ATF6" = state_palette["ATF6-high"])) +
  labs(title = "Mean UPR Arm Scores by Reostat State", 
       subtitle = "Colorful bars show relative arm activation levels",
       x = "UPR Reostat State", y = "Mean Activation Score", fill = "UPR Arm") +
  theme_publication(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    axis.text.y = element_text(size = 11),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 11, color = "grey30"),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 11),
    legend.position = "bottom"
  ) +
  coord_cartesian(ylim = c(0, NA)) +
  guides(fill = guide_legend(nrow = 1))

# Enhanced Figure 3 layout
fig3 <- (p3A | p3B) / (p3C | p3D) + 
  plot_annotation(
    tag_levels = 'A',
    title = "Figure 3: UPR-Stress Integration and Cellular Dynamics",
    subtitle = "Coordination of UPR pathways with stress responses and temporal dynamics",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "grey30")
    )
  ) +
  plot_layout(heights = c(1, 1))

# Save as PDF
ggsave(out_fig3, fig3, width = 16, height = 12, units = "in", dpi = 300)
log_message("Figure 3 PDF saved:", out_fig3)

# Save as high-resolution PNG
png_fig3 <- gsub("\\.pdf$", "_1200dpi.png", out_fig3)
ggsave(png_fig3, fig3, width = 16, height = 12, units = "in", dpi = 1200, bg = "white")
log_message("Figure 3 PNG saved:", png_fig3)

# Figure summary statistics
log_message("=== FIGURE SUMMARY ===")
log_message("Total cell count:", nrow(plot_df))
log_message("UPR state distribution:")
upr_dist <- table(plot_df$UPR_State)
for(i in 1:length(upr_dist)) {
  log_message(paste0("  ", names(upr_dist)[i], ": ", upr_dist[i], " (", 
                     round(upr_dist[i]/sum(upr_dist)*100, 1), "%)"))
}
log_message("Correlations (Spearman):")
log_message(paste0("  PERK-IRE1: r=", round(cor_perk_ire1$cor, 3), " ", format_p_value(cor_perk_ire1$p)))
log_message(paste0("  PERK-ATF6: r=", round(cor_perk_atf6$cor, 3), " ", format_p_value(cor_perk_atf6$p)))
log_message(paste0("  IRE1-ATF6: r=", round(cor_ire1_atf6$cor, 3), " ", format_p_value(cor_ire1_atf6$p)))

log_message("All figures successfully generated. Nature-quality multi-panel figures ready!")


