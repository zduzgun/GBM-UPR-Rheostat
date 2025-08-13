# ==============================================================================
# 07_10_rheostat_validation_on_trajectory.R (v2 - Corrected Input File)
#
# Purpose: To directly test the UPR Rheostat hypothesis on the Neftel data
#          by visualizing UPR dynamics along the malignant cell trajectory.
#          FIX: Using the correct input RDS file from step 07_07.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(monocle3)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(patchwork) 
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
# ================================ DÜZELTME BURADA =================================
# Doðru girdi dosyasý, 07_07'den gelen CDS nesnesidir.
input_cds_file <- "07_07_malignant_monocle_cds.rds"
# ===============================================================================
output_plot_file <- "results/figures/07_10_rheostat_validation_plots.pdf"
output_data_file <- "results/tables/07_10_rheostat_summary_by_pseudotime.csv"


# --- 3. Load Processed CDS Object ---
cat("Step 3: Loading the CDS object with pseudotime and scores...\n")
if (!file.exists(input_cds_file)) {
  stop("Input file not found: ", input_cds_file, ". Please run 07_07 script first.")
}
cds <- readRDS(input_cds_file)
cat("CDS object loaded successfully.\n")


# --- 4. Prepare Data for Plotting ---
cat("Step 4: Preparing data for visualization...\n")
cell_data <- as.data.frame(colData(cds))

required_cols <- c("pseudotime", "IRE1_score", "PERK_score", "ATF6_score")
if (!all(required_cols %in% colnames(cell_data))) {
  missing <- required_cols[!required_cols %in% colnames(cell_data)]
  stop("Required columns are missing from the CDS metadata: ", paste(missing, collapse=", "))
}

if (!"UPR_Reostat_SD" %in% colnames(cell_data)) {
  cat("  - Calculating UPR_Reostat_SD on the fly...\n")
  upr_matrix <- as.matrix(cell_data[, c("IRE1_score", "PERK_score", "ATF6_score")])
  cell_data$UPR_Reostat_SD <- apply(upr_matrix, 1, sd, na.rm = TRUE)
}

data_long <- cell_data %>%
  select(pseudotime, IRE1_score, PERK_score, ATF6_score) %>%
  pivot_longer(cols = c(IRE1_score, PERK_score, ATF6_score),
               names_to = "upr_arm",
               values_to = "score")


# --- 5. Create Plots ---
cat("Step 5: Generating plots to test the rheostat hypothesis...\n")

p1 <- ggplot(data_long, aes(x = pseudotime, y = score, color = upr_arm)) +
  geom_smooth(method = "loess", se = TRUE, span = 0.5) +
  scale_color_manual(values = c("IRE1_score" = "#E31A1C", "PERK_score" = "#1F78B4", "ATF6_score" = "#33A02C")) +
  labs(
    title = "Dynamics of UPR Arms Along Malignant Trajectory",
    subtitle = "Do the arms show different patterns over pseudotime?",
    x = "Pseudotime (from Stem-like to Differentiated)",
    y = "UPR Arm Activity Score",
    color = "UPR Arm"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

p2 <- ggplot(cell_data, aes(x = pseudotime, y = UPR_Reostat_SD)) +
  geom_point(alpha = 0.1, color = "grey50") +
  geom_smooth(method = "loess", color = "purple", se = TRUE, span = 0.5) +
  labs(
    title = "UPR Rheostat Index (Variability) vs. Pseudotime",
    subtitle = "Does the imbalance between UPR arms change during progression?",
    x = "Pseudotime",
    y = "UPR Rheostat Index (StDev)"
  ) +
  theme_minimal()


# --- 6. Save Combined Plot ---
cat("Step 6: Saving combined plots to PDF...\n")
combined_plot <- p1 / p2

pdf(output_plot_file, width = 10, height = 12)
print(combined_plot)
dev.off()


# --- 7. Save Summary Data (Optional) ---
cat("Step 7: Saving summary data...\n")
summary_data <- cell_data %>%
  mutate(pseudotime_bin = cut_number(pseudotime, n = 20)) %>%
  group_by(pseudotime_bin) %>%
  summarise(
    mean_pseudotime = mean(pseudotime),
    mean_IRE1 = mean(IRE1_score),
    mean_PERK = mean(PERK_score),
    mean_ATF6 = mean(ATF6_score),
    mean_Reostat_SD = mean(UPR_Reostat_SD)
  )
write.csv(summary_data, file = output_data_file, row.names = FALSE)


cat("Analysis complete. Check for the output PDF:", output_plot_file, "\n")