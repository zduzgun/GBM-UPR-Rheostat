#!/usr/bin/env Rscript

# Figure 3: UPR-Stress Integration and Cellular Dynamics Generator
# Standalone script to generate Figure 3 from glioblastoma_with_corrected_rheostat.rds

options(stringsAsFactors = FALSE, scipen = 999)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
  library(scales)
  library(viridis)
  library(corrplot)
})

# Utility functions
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

# Enhanced publication theme
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

# Helper functions for object handling
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
    for (nm in c("umap", "UMAP", "tsne", "pca")) {
      if (nm %in% names(x@reductions)) { red_name <- nm; break }
    }
    if (is.null(red_name)) {
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

# UPR state classification function
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

# ----------------------
# MAIN EXECUTION
# ----------------------

# Input and output files
input_rds <- "glioblastoma_with_corrected_rheostat.rds"
out_fig3 <- "Figure3_UPR_Stress_Integration.pdf"

# Check if input file exists
if (!file.exists(input_rds)) {
  stop_with_hint(sprintf("Expected RDS file not found: %s", input_rds))
}

log_message("Loading object:", input_rds)
obj <- readRDS(input_rds)

# Get metadata and coordinates
meta <- get_metadata_df(obj)
emb_list <- get_embeddings_df(obj)
emb <- emb_list$emb
obj <- emb_list$obj

log_message("Available columns in metadata:")
log_message(paste(colnames(meta), collapse = ", "))

# Identify UPR and other important columns
cluster_col <- choose_first_available(colnames(meta), c("^seurat_clusters$", "cluster$", "cluster_id$", "snn_res", "leiden", "louvain"))
if (is.na(cluster_col)) cluster_col <- choose_first_available(colnames(meta), c("cluster", "group", "subcluster"))
if (is.na(cluster_col)) cluster_col <- "cluster"

perk_col <- choose_first_available(colnames(meta), c("perk", "eif2ak3", "perk_score", "UPR_PERK"))
ire1_col <- choose_first_available(colnames(meta), c("ire1", "ern1", "xbp1", "ire1_score", "UPR_IRE1"))
atf6_col <- choose_first_available(colnames(meta), c("atf6", "atf6_score", "UPR_ATF6"))
reostat_col <- choose_first_available(colnames(meta), c("reostat", "upr_reostat", "upr_reostat_score", "UPR_Reostat"))

log_message("Found columns - PERK:", perk_col, "IRE1:", ire1_col, "ATF6:", atf6_col, "Reostat:", reostat_col)

# Ensure UPR columns exist
if (is.na(perk_col) || is.na(ire1_col) || is.na(atf6_col)) {
  stop_with_hint("UPR arms (PERK, IRE1, ATF6) not found in metadata.")
}

# Standardize UPR columns
meta <- meta %>% mutate(
  UPR_PERK = .data[[perk_col]],
  UPR_IRE1 = .data[[ire1_col]],
  UPR_ATF6 = .data[[atf6_col]]
)

# Add reostat if available
if (!is.na(reostat_col)) {
  meta$UPR_Reostat <- meta[[reostat_col]]
}

# Classify UPR states
meta <- meta %>% mutate(
  UPR_State = classify_reostat_state(., perk = "UPR_PERK", ire1 = "UPR_IRE1", atf6 = "UPR_ATF6")
)

# Merge with coordinates
plot_df <- emb %>% left_join(meta, by = "cell_id")

# Color palettes
state_palette <- c(
  "All-low" = "#CCCCCC",
  "PERK-high" = "#E31A1C", 
  "IRE1-high" = "#1F78B4",
  "ATF6-high" = "#33A02C",
  "Multi-arm-high" = "#9A3192",
  "Intermediate/Mixed" = "#FF7F00"
)

# Enhanced stress module detection
stress_like_cols <- grep("hypoxia|oxidative|stress|emt|apoptosis|dna_damage|unfolded|proteasome|inflammation|nfkb|tnfa|ifn|il6|glycolysis|oxidative_phosphorylation|cholesterol|lipid|autophagy|angiogenesis|stemness|differentiation",
                         colnames(meta), ignore.case = TRUE, value = TRUE)

sel_cols <- unique(c("UPR_PERK", "UPR_IRE1", "UPR_ATF6", intersect(stress_like_cols, colnames(meta))))

log_message("Stress-related columns found:", paste(intersect(stress_like_cols, colnames(meta)), collapse = ", "))

# ----------------------
# Panel A: Correlation matrix
# ----------------------
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
    select(where(~ !all(is.na(.x))))
  
  if (ncol(sub) < 2) {
    return(ggplot() + theme_void() + labs(title = "Insufficient valid data for correlation"))
  }
  
  cm <- suppressWarnings(cor(sub, use = "pairwise.complete.obs", method = "pearson"))
  
  # Calculate p-values safely
  n_vars <- ncol(sub)
  p_matrix <- matrix(1, n_vars, n_vars)
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

# ----------------------
# Panel B: UPR arm balance ratios
# ----------------------
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
  filter(is.finite(Value) & Value > 0)

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

# ----------------------
# Panel C: UPR reostat manifold
# ----------------------
color_col <- if (!is.na(reostat_col)) "UPR_Reostat" else "UPR_State"

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

# Panel D removed as requested

# ----------------------
# Combine panels into Figure 3 (3 panels: A, B, C)
# ----------------------
fig3 <- (p3A | p3B) / (p3C) + 
  plot_annotation(
    tag_levels = list(c('a', 'b', 'c')),
    title = "Figure 3: UPR-Stress Integration and Cellular Dynamics",
    subtitle = "Coordination of UPR pathways with stress responses and temporal dynamics",
    theme = theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "grey30"),
      plot.tag = element_text(size = 20, face = "bold")
    )
  ) +
  plot_layout(heights = c(1, 1))

# Save Figure 3
ggsave(out_fig3, fig3, width = 16, height = 12, units = "in", dpi = 300)
log_message("Figure 3 PDF saved:", out_fig3)

# Save as high-resolution PNG
png_fig3 <- gsub("\\.pdf$", "_1200dpi.png", out_fig3)
ggsave(png_fig3, fig3, width = 16, height = 12, units = "in", dpi = 1200, bg = "white")
log_message("Figure 3 PNG saved:", png_fig3)

# ----------------------
# Export data tables as CSV
# ----------------------

# 1. UPR correlation matrix data
if (length(sel_cols) >= 3) {
  cor_data <- meta %>% 
    select(all_of(sel_cols)) %>% 
    mutate(across(everything(), as.numeric)) %>%
    select(where(~ !all(is.na(.x))))
  
  if (ncol(cor_data) >= 2) {
    cm <- suppressWarnings(cor(cor_data, use = "pairwise.complete.obs", method = "pearson"))
    write.csv(cm, "Figure3A_correlation_matrix.csv", row.names = TRUE)
    log_message("Correlation matrix saved: Figure3A_correlation_matrix.csv")
  }
}

# 2. UPR balance ratios data
write.csv(upr_balance, "Figure3B_upr_balance_ratios.csv", row.names = FALSE)
log_message("UPR balance ratios saved: Figure3B_upr_balance_ratios.csv")

# 3. Plot coordinates with UPR data
coordinates_data <- plot_df %>%
  select(cell_id, Dim1, Dim2, UPR_PERK, UPR_IRE1, UPR_ATF6, UPR_State, 
         all_of(cluster_col), any_of("UPR_Reostat"))
write.csv(coordinates_data, "Figure3C_coordinates_upr_data.csv", row.names = FALSE)
log_message("Coordinates and UPR data saved: Figure3C_coordinates_upr_data.csv")
write.csv(meta, "Figure3_complete_metadata.csv", row.names = FALSE)
log_message("Complete metadata saved: Figure3_complete_metadata.csv")

# ----------------------
# Summary statistics
# ----------------------
log_message("=== FIGURE 3 SUMMARY ===")
log_message("Total cell count:", nrow(plot_df))
log_message("UPR state distribution:")
upr_dist <- table(plot_df$UPR_State)
for(i in 1:length(upr_dist)) {
  log_message(paste0("  ", names(upr_dist)[i], ": ", upr_dist[i], " (", 
                     round(upr_dist[i]/sum(upr_dist)*100, 1), "%)"))
}

log_message("Available stress-related columns:", paste(intersect(stress_like_cols, colnames(meta)), collapse = ", "))
log_message("Figure 3 generation completed successfully!")