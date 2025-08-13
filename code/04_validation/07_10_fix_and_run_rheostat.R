# ==============================================================================
# 07_10_fix_and_run_rheostat.R (v4 - Fixed Directory Issue)
#
# Purpose: 1. Load the trajectory CDS object (missing scores) and the original
#             scored Seurat object.
#          2. Merge the UPR scores into the CDS object by matching cell barcodes.
#          3. Directly test the UPR Rheostat hypothesis by visualizing the
#             dynamics of UPR arms along the malignant cell trajectory.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)       # Seurat nesnesini okumak için eklendi
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork)
  library(tibble)     # rownames_to_column için eklendi
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
cds_trajectory_file <- "07_07_malignant_monocle_cds.rds"
seurat_with_scores_file <- "07_03_neftel_with_scores.rds" # <-- Teyit ettiðiniz dosya adý
output_plot_file <- "07_10_fix_and_run_rheostat_results/figures/07_10_rheostat_validation_plots.pdf"
output_data_file <- "07_10_fix_and_run_rheostat_results/tables/07_10_rheostat_summary_by_pseudotime.csv"

# --- 2.1. Create Output Directories if They Don't Exist ---
cat("Step 2.1: Creating output directories...\n")
dir.create("07_10_fix_and_run_rheostat_results", showWarnings = FALSE, recursive = TRUE)
dir.create("07_10_fix_and_run_rheostat_results/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("07_10_fix_and_run_rheostat_results/tables", showWarnings = FALSE, recursive = TRUE)
cat("  - Output directories created/verified.\n")

# --- 3. Load Both Data Objects ---
cat("Step 3: Loading data objects...\n")
if (!file.exists(cds_trajectory_file)) stop("Trajectory CDS file not found: ", cds_trajectory_file)
if (!file.exists(seurat_with_scores_file)) stop("Seurat scores file not found: ", seurat_with_scores_file)

cds <- readRDS(cds_trajectory_file)
seurat_obj <- readRDS(seurat_with_scores_file)
cat("  - Trajectory CDS loaded:", ncol(cds), "cells.\n")
cat("  - Scored Seurat object loaded:", ncol(seurat_obj), "cells.\n")

# --- 4. FIX: Merge Missing UPR Scores into the CDS Object ---
cat("Step 4: Merging missing UPR scores into the CDS object...\n")

# Metadata'larý (üst verileri) al
cds_metadata <- as.data.frame(colData(cds)) %>% rownames_to_column("cell_id")
seurat_metadata <- seurat_obj@meta.data %>% rownames_to_column("cell_id")

# Seurat nesnesindeki skor sütunlarýný bul (hem _score hem de _score1 formatýný ara)
score_cols <- grep("IRE1|PERK|ATF6", colnames(seurat_metadata), value = TRUE)
if (length(score_cols) == 0) stop("No UPR score columns found in the Seurat object!")
cat("  - Found score columns in Seurat object:", paste(score_cols, collapse=", "), "\n")

# Sadece skor sütunlarýný ve hücre ID'sini seç
seurat_scores <- seurat_metadata %>% select(cell_id, all_of(score_cols))

# CDS metadata'sýna skorlarý hücre ID'sine göre ekle (sol birleþtirme)
merged_metadata <- left_join(cds_metadata, seurat_scores, by = "cell_id")

# Hata kontrolü: birleþme sonrasý satýr sayýsýnýn deðiþmediðinden emin ol
if (nrow(merged_metadata) != nrow(cds_metadata)) {
  stop("Metadata merge failed! Row count changed.")
}

# Birleþtirilmiþ ve düzeltilmiþ metadata'yý CDS nesnesine geri koy
# Önce satýr isimlerini (hücre ID'leri) geri ayarla
rownames(merged_metadata) <- merged_metadata$cell_id
colData(cds) <- DataFrame(merged_metadata, row.names = merged_metadata$cell_id)

cat("  - Scores successfully merged into the CDS object.\n")

# --- 5. Prepare Data for Plotting (Using Corrected Names) ---
cat("Step 5: Preparing data for visualization...\n")
cell_data <- as.data.frame(colData(cds)) # Artýk skorlarý içeren güncel veri

# Ýsimleri standart hale getir (_score1 varsa _score yap)
colnames(cell_data) <- sub("_score1$", "_score", colnames(cell_data))
colnames(cell_data) <- sub("_Score1$", "_score", colnames(cell_data))

# Gerekli sütunlarýn varlýðýna son kez kontrol et
required_cols <- c("pseudotime", "IRE1_score", "PERK_score", "ATF6_score")
if (!all(required_cols %in% colnames(cell_data))) {
  stop("Even after merge, required columns are missing. Check column names.")
}

if (!"UPR_Reostat_SD" %in% colnames(cell_data)) {
  upr_matrix <- as.matrix(cell_data[, required_cols[2:4]])
  cell_data$UPR_Reostat_SD <- apply(upr_matrix, 1, sd, na.rm = TRUE)
}

data_long <- cell_data %>%
  select(pseudotime, all_of(required_cols[2:4])) %>%
  pivot_longer(cols = -pseudotime, names_to = "upr_arm", values_to = "score")

# --- 6. Create Plots ---
cat("Step 6: Generating plots to test the rheostat hypothesis...\n")
# (Bu kýsým öncekiyle ayný, çünkü artýk veri doðru)
p1 <- ggplot(data_long, aes(x = pseudotime, y = score, color = upr_arm)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.5) +
  scale_color_manual(values = c("IRE1_score" = "#E31A1C", "PERK_score" = "#1F78B4", "ATF6_score" = "#33A02C")) +
  labs(title = "Dynamics of UPR Arms Along Malignant Trajectory", x = "Pseudotime", y = "UPR Arm Activity Score") +
  theme_minimal() + theme(legend.position = "bottom")

p2 <- ggplot(cell_data, aes(x = pseudotime, y = UPR_Reostat_SD)) +
  geom_point(alpha = 0.1, color = "grey50") +
  geom_smooth(method = "loess", color = "purple", se = TRUE, span = 0.5) +
  labs(title = "UPR Rheostat Index (Variability) vs. Pseudotime", x = "Pseudotime", y = "UPR Rheostat Index (StDev)") +
  theme_minimal()

# --- 7. Save Combined Plot and Data ---
cat("Step 7: Saving outputs...\n")
combined_plot <- p1 / p2
pdf(output_plot_file, width = 10, height = 12)
print(combined_plot)
dev.off()

summary_data <- cell_data %>%
  mutate(pseudotime_bin = cut_number(pseudotime, n = 20)) %>%
  group_by(pseudotime_bin) %>%
  summarise(across(c(pseudotime, IRE1_score, PERK_score, ATF6_score, UPR_Reostat_SD), mean, na.rm=TRUE))
write.csv(summary_data, file = output_data_file, row.names = FALSE)

cat("\nAnalysis complete. Check for the output PDF:", output_plot_file, "\n")