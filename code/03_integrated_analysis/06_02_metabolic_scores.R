#################################################################
# SCRIPT 06_02: METABOLIC PATHWAY SCORE CALCULATION
# Purpose: To calculate module scores for key metabolic pathways
#          (Glycolysis and OXPHOS) for each cell.
# Input:   glioblastoma_with_all_scores.rds
# Output:  glioblastoma_with_final_scores.rds
#################################################################

# Load required libraries
library(Seurat)
library(dplyr)

# --- Load Data ---
# Load the Seurat object containing UPR and Functional Axis scores
print("Reading Seurat object: glioblastoma_with_all_scores.rds")
seurat_obj <- readRDS("glioblastoma_with_all_scores.rds")
print("Data loading complete.")

# --- DEFINE GENE LISTS ---
# Metabolic Pathway Gene Sets (from MSigDB Hallmark sets)
glycolysis_genes <- c('ALDOA', 'BPGM', 'ENO1', 'ENO2', 'GAPDH', 'GPI', 'HK1', 'HK2', 'LDHA', 'PFKL', 'PFKM', 'PFKP', 'PGAM1', 'PGK1', 'PKM', 'TPI1')
oxphos_genes <- c('ATP5F1A', 'ATP5F1B', 'ATP5MC1', 'ATP5MC2', 'ATP5MC3', 'COX4I1', 'COX5A', 'COX5B', 'COX6C', 'COX7B', 'NDUFA1', 'NDUFB1', 'NDUFC2', 'SDHA', 'SDHB', 'UQCRB', 'UQCRC1', 'UQCRQ')

# --- Filter gene lists to keep only genes present in the Seurat object ---
genes_in_object <- rownames(seurat_obj)

filter_genes <- function(gene_list, all_genes) {
  return(list(intersect(gene_list, all_genes)))
}

glycolysis_genes_filtered <- filter_genes(glycolysis_genes, genes_in_object)
oxphos_genes_filtered <- filter_genes(oxphos_genes, genes_in_object)

# Print informative messages
print("Metabolic gene lists have been filtered. Number of genes found:")
print(paste("Glycolysis:", length(glycolysis_genes_filtered[[1]]), "/", length(glycolysis_genes)))
print(paste("OXPHOS:", length(oxphos_genes_filtered[[1]]), "/", length(oxphos_genes)))


# --- CALCULATE MODULE SCORES ---
print("Calculating Metabolic Pathway scores...")
seurat_obj <- AddModuleScore(seurat_obj, features = glycolysis_genes_filtered, name = "Glycolysis_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = oxphos_genes_filtered, name = "OXPHOS_Score")
print("Scoring is complete.")


# --- SAVE THE UPDATED DATA ---
# This object now contains all scores (UPR, Functional, Metabolic)
print("Saving the updated Seurat object: glioblastoma_with_final_scores.rds")
saveRDS(seurat_obj, "glioblastoma_with_final_scores.rds")

print("Process completed successfully.")