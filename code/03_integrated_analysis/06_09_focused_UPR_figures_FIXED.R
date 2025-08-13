#!/usr/bin/env Rscript

# 06_09_focused_UPR_figures_FIXED.R
# Focused figures for UPR Reostat hypothesis in GBM heterogeneity
# CORRECTED to fix scientific, layout, and coding errors.

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ComplexHeatmap)
  library(patchwork)
  library(RColorBrewer)
  library(viridis)
  library(circlize) # DÜZELTME: circlize::colorRamp2 için açýkça eklendi
})

message("=== Focused UPR Reostat Analysis (Corrected) ===")
message("Core research: GBM heterogeneity and UPR dynamic regulation")

# Load data
# DÜZELTME: Veriyi bu script içinde düzelteceðimiz için orijinal dosyayý yüklüyoruz.
# 'glioblastoma_with_corrected_reostat.rds' yerine 'glioblastoma_with_final_scores.rds' kullanýlmasý gerekebilir.
# Bizde bu dosya olmadýðý için, kodun çalýþmasý adýna mevcut olaný yüklüyoruz.
tryCatch({
    glioblastoma <- readRDS("glioblastoma_with_final_scores.rds")
}, error = function(e) {
    message("Could not find 'glioblastoma_with_final_scores.rds', trying 'glioblastoma_with_corrected_reostat.rds'")
    glioblastoma <- readRDS("glioblastoma_with_corrected_reostat.rds")
})


# Verify UPR scores are available
upr_scores <- c("IRE1_score", "PERK_score", "ATF6_score")
if(!all(upr_scores %in% colnames(glioblastoma@meta.data))) {
  stop("UPR scores not found in metadata.")
}

# Professional theme
theme_pub <- theme_minimal() +
  theme(
    text = element_text(size = 11),
    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 9),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.minor = element_blank(),
    strip.text = element_text(size = 10, face = "bold")
  )

# Color palette for clusters
cluster_colors <- RColorBrewer::brewer.pal(9, "Set1")
names(cluster_colors) <- 0:8

# === FIGURE 1: GBM Heterogeneity Overview ===
message("Creating Figure 1: GBM Intratumoral Heterogeneity...")

pdf("Figure1_GBM_Heterogeneity.pdf", width = 14, height = 10)

p1a <- DimPlot(glioblastoma, group.by = "seurat_clusters",
               cols = cluster_colors, pt.size = 0.4, label = TRUE,
               label.size = 4, repel = TRUE) +
  labs(title = "Intratumoral Cellular Heterogeneity",
       subtitle = "9 distinct cell states in GBM ecosystem") +
  theme_pub +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank()) +
  guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))

phenotype_scores <- c("Proliferation_Score1", "Mesenchymal_Score1", "Hypoxia_Score1")
phenotype_names <- c("Proliferative", "Mesenchymal", "Hypoxic")

phenotype_plots <- list()
for(i in seq_along(phenotype_scores)) {
  phenotype_plots[[i]] <- FeaturePlot(glioblastoma, features = phenotype_scores[i],
                                     pt.size = 0.3, cols = c("lightgray", "darkred")) +
    labs(title = paste(phenotype_names[i], "Program")) +
    theme_pub +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(), plot.title = element_text(size = 11))
}

cluster_stats <- glioblastoma@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(cell_count = n(), .groups = 'drop')

p1c <- ggplot(cluster_stats, aes(x = factor(seurat_clusters), y = cell_count,
                                fill = factor(seurat_clusters))) +
  geom_col(alpha = 0.8) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "Cellular Composition",
       subtitle = "Cell count per cluster",
       x = "Cell Cluster", y = "Number of Cells") +
  theme_pub +
  guides(fill = "none")

## DÜZELTME: Figür 1 düzeni, PDF'deki görünüme uyacak þekilde basitleþtirildi ve düzeltildi.
right_panel <- (phenotype_plots[[1]] / phenotype_plots[[2]] / phenotype_plots[[3]])
combined_fig1 <- (p1a | right_panel) / p1c +
  plot_layout(heights = c(2.5, 1)) +
  plot_annotation(
    title = "Figure 1: Glioblastoma Intratumoral Heterogeneity",
    subtitle = "Single-cell analysis reveals diverse cellular phenotypes and states",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

print(combined_fig1)
dev.off()

# === FIGURE 2: UPR Reostat Hypothesis Testing ===
message("Creating Figure 2: UPR Reostat Model...")

pdf("Figure2_UPR_Reostat_Hypothesis.pdf", width = 16, height = 12)

upr_data <- glioblastoma@meta.data %>%
  select(seurat_clusters, IRE1_score, PERK_score, ATF6_score) %>%
  reshape2::melt(id.vars = "seurat_clusters", variable.name = "UPR_Pathway",
                 value.name = "Activity_Score")

p2a <- ggplot(upr_data, aes(x = factor(seurat_clusters), y = Activity_Score,
                           fill = UPR_Pathway)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.15,
               outlier.size = 0.5, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", position = position_dodge(0.8),
               size = 1.5, color = "black") +
  scale_fill_manual(values = c("IRE1_score" = "#E31A1C",
                              "PERK_score" = "#1F78B4",
                              "ATF6_score" = "#33A02C"),
                   labels = c("IRE1?", "PERK", "ATF6"),
                   name = "UPR Pathway") +
  labs(title = "UPR Pathway Activities Across Cell States",
       subtitle = "Evidence for differential pathway activation patterns",
       x = "Cell Cluster", y = "UPR Activity Score") +
  theme_pub +
  theme(legend.position = "bottom")

cor_matrix <- cor(glioblastoma@meta.data[, upr_scores], use = "complete.obs")
cor_df <- expand.grid(Pathway1 = rownames(cor_matrix), Pathway2 = colnames(cor_matrix))
cor_df$Correlation <- as.vector(cor_matrix)
cor_df$Pathway1 <- factor(cor_df$Pathway1, levels = rev(upr_scores))

p2b <- ggplot(cor_df, aes(Pathway2, Pathway1, fill = Correlation)) +
  ## DÜZELTME: Terminal uyarýsýný gidermek için 'size' yerine 'linewidth' kullanýldý.
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = sprintf("r = %.2f", Correlation)),
            size = 4, fontface = "bold") +
  scale_fill_gradient2(low = "#1F78B4", mid = "white", high = "#E31A1C",
                      midpoint = 0, name = "Correlation\nCoefficient",
                      limits = c(-1, 1)) +
  scale_x_discrete(labels = c("IRE1?", "PERK", "ATF6")) +
  scale_y_discrete(labels = c("ATF6", "PERK", "IRE1?")) +
  labs(title = "UPR Pathway Correlations",
       subtitle = "Moderate correlations support reostat model",
       x = "", y = "") +
  theme_pub +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## DÜZELTME: KRÝTÝK BÝLÝMSEL HATA DÜZELTÝLDÝ
# Reostat indeksi, hatalý olan Varyasyon Katsayýsý (CV) yerine,
# bilimsel olarak daha saðlam olan Standart Sapma (SD) olarak hesaplandý.
upr_matrix <- as.matrix(glioblastoma@meta.data[, upr_scores])
reostat_sd <- apply(upr_matrix, 1, function(x) sd(x, na.rm = TRUE))
glioblastoma$UPR_Reostat_SD <- reostat_sd

# Düzeltilmiþ indeksi gösteren FeaturePlot
p2c <- FeaturePlot(glioblastoma, features = "UPR_Reostat_SD",
                   pt.size = 0.3, cols = c("lightblue", "navy", "red")) +
  labs(title = "UPR Reostat Index (Standard Deviation)",
       subtitle = "Pathway variability within cells") +
  theme_pub +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.grid = element_blank())

# Düzeltilmiþ indeksi kümelere göre özetleyen bar plot
reostat_by_cluster <- glioblastoma@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    mean_sd = mean(UPR_Reostat_SD, na.rm = TRUE),
    se_sd = sd(UPR_Reostat_SD, na.rm = TRUE) / sqrt(n()),
    .groups = 'drop'
  )

p2d <- ggplot(reostat_by_cluster, aes(x = factor(seurat_clusters), y = mean_sd,
                                     fill = factor(seurat_clusters))) +
  geom_col(alpha = 0.8, width = 0.7) +
  geom_errorbar(aes(ymin = mean_sd - se_sd, ymax = mean_sd + se_sd),
                width = 0.3, alpha = 0.8) +
  scale_fill_manual(values = cluster_colors) +
  labs(title = "UPR Reostat by Cell State",
       subtitle = "Mean ± SE Standard Deviation",
       x = "Cell Cluster", y = "UPR Reostat Index (SD)") +
  theme_pub +
  guides(fill = "none")


layout_fig2 <- "
AAAABBBB
AAAABBBB
CCCCDDDD
"
combined_fig2 <- p2a + p2b + p2c + p2d +
  plot_layout(design = layout_fig2) +
  plot_annotation(
    title = "Figure 2: UPR Reostat Hypothesis Validation",
    subtitle = "Dynamic UPR pathway regulation supports flexible stress response model",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 12, hjust = 0.5))
  )

print(combined_fig2)
dev.off()


# === FIGURE 3: UPR-Stress Relationship ===
message("Creating Figure 3: UPR-Stress Integration...")

## DÜZELTME: Boþ PDF dosyasý oluþturan gereksiz pdf() çaðrýsý kaldýrýldý.
# pdf("Figure3_UPR_Stress_Integration.pdf", width = 14, height = 10)

p3a <- ggplot(glioblastoma@meta.data,
              aes(x = Hypoxia_Score1, y = IRE1_score, color = factor(seurat_clusters))) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  labs(title = "IRE1? vs Hypoxic Stress",
       x = "Hypoxia Score", y = "IRE1? Activity") +
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 3)))

p3b <- ggplot(glioblastoma@meta.data,
              aes(x = Proliferation_Score1, y = PERK_score, color = factor(seurat_clusters))) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  labs(title = "PERK vs Proliferative State",
       x = "Proliferation Score", y = "PERK Activity") +
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 3)))

p3c <- ggplot(glioblastoma@meta.data,
              aes(x = Mesenchymal_Score1, y = ATF6_score, color = factor(seurat_clusters))) +
  geom_point(alpha = 0.6, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
  scale_color_manual(values = cluster_colors, name = "Cluster") +
  labs(title = "ATF6 vs Mesenchymal Program",
       x = "Mesenchymal Score", y = "ATF6 Activity") +
  theme_pub +
  guides(color = guide_legend(override.aes = list(size = 3)))

all_scores <- c("IRE1_score", "PERK_score", "ATF6_score",
                "Hypoxia_Score1", "Proliferation_Score1", "Mesenchymal_Score1")
score_labels <- c("IRE1?", "PERK", "ATF6", "Hypoxia", "Proliferation", "Mesenchymal")

cor_all <- cor(glioblastoma@meta.data[, all_scores], use = "complete.obs")
rownames(cor_all) <- colnames(cor_all) <- score_labels

col_fun <- circlize::colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red"))

p3d_heatmap <- ComplexHeatmap::Heatmap(
  cor_all,
  name = "Correlation",
  col = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if(i != j) {
      grid::grid.text(sprintf("%.2f", cor_all[i, j]), x, y,
                     gp = grid::gpar(fontsize = 9))
    }
  },
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  column_title = "UPR-Phenotype Correlation Matrix",
  column_title_gp = grid::gpar(fontsize = 12, fontface = "bold")
)

combined_scatter <- p3a / p3b / p3c

pdf("Figure3_UPR_Stress_Integration_combined.pdf", width = 16, height = 12)
print(combined_scatter)
dev.off()

pdf("Figure3D_Correlation_Matrix.pdf", width = 8, height = 8)
draw(p3d_heatmap)
dev.off()

## DÜZELTME: Fazladan ve hataya neden olan dev.off() çaðrýsý kaldýrýldý.
# dev.off()

# === SUMMARY STATISTICS ===
message("Calculating summary statistics...")

upr_correlations <- cor(glioblastoma@meta.data[, upr_scores], use = "complete.obs")
message("UPR Pathway Correlations:")
message(paste("IRE1-PERK:", round(upr_correlations["IRE1_score", "PERK_score"], 3)))
message(paste("IRE1-ATF6:", round(upr_correlations["IRE1_score", "ATF6_score"], 3)))
message(paste("PERK-ATF6:", round(upr_correlations["PERK_score", "ATF6_score"], 3)))

## DÜZELTME: Özet istatistikler, düzeltilmiþ Reostat Ýndeksi (SD) kullanacak þekilde güncellendi.
cluster_upr_summary <- glioblastoma@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    IRE1_mean = mean(IRE1_score, na.rm = TRUE),
    PERK_mean = mean(PERK_score, na.rm = TRUE),
    ATF6_mean = mean(ATF6_score, na.rm = TRUE),
    Reostat_SD_mean = mean(UPR_Reostat_SD, na.rm = TRUE),
    .groups = 'drop'
  )

write.csv(cluster_upr_summary, "UPR_cluster_summary_corrected.csv", row.names = FALSE)

message("\n=== FOCUSED UPR REOSTAT FIGURES COMPLETED (CORRECTED) ===")
message("Generated files:")
message("• Figure1_GBM_Heterogeneity.pdf - Intratumoral heterogeneity")
message("• Figure2_UPR_Reostat_Hypothesis.pdf - Core hypothesis testing")
message("• Figure3_UPR_Stress_Integration_combined.pdf - UPR-phenotype relationships")
message("• Figure3D_Correlation_Matrix.pdf - Comprehensive correlation analysis")
message("• UPR_cluster_summary_corrected.csv - Corrected statistical summary")
message("\nAll figures and calculations have been corrected for scientific accuracy and technical errors.")