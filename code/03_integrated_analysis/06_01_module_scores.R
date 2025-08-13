#################################################################
# SCRIPT 06_01: MODULE SCORE CALCULATION
# Purpose: To calculate module scores for UPR arms and key
#          functional axes for each cell.
# Input:   ../glioblastoma_processed.rds
# Output:  glioblastoma_with_all_scores.rds
#################################################################

# Load required libraries
library(Seurat)
library(dplyr)

# --- Load Data ---
print("Reading Seurat object: ../glioblastoma_processed.rds")
seurat_obj <- readRDS("../glioblastoma_processed.rds")
print("Data loading complete.")

# --- DEFINE GENE LISTS ---

# UPR Pathway Gene Sets
ire1_genes <- c('XBP1', 'DNAJB9', 'EDEM1', 'HSPA5', 'PDIA4', 'ERO1A', 'ERN1')
perk_genes <- c('ATF4', 'DDIT3', 'PPP1R15A', 'ASNS', 'TRIB3', 'EIF2AK3')
atf6_genes <- c('ATF6', 'HSPA5', 'HSP90B1', 'CALR')

# Functional Axis Gene Sets
proliferation_genes <- c('MKI67', 'TOP2A', 'CDK1', 'KIF15', 'H1-3')
mesenchymal_genes <- c('SPP1', 'CHI3L1', 'FABP5', 'POSTN', 'VIM', 'FN1')
stemness_genes <- c('MYCN', 'MAGEA3', 'MAGEA12', 'SOX2', 'PROM1') # PROM1 = CD133
hypoxia_genes <- c('VEGFA', 'NDRG1', 'HILPDA', 'ADM', 'ANKRD37')

# --- Filter gene lists to keep only genes present in the Seurat object ---

# Get a list of all genes present in the Seurat object
genes_in_object <- rownames(seurat_obj)

# Function to filter a gene list against the genes in the object
filter_genes <- function(gene_list, all_genes) {
  return(list(intersect(gene_list, all_genes)))
}

ire1_genes_filtered <- filter_genes(ire1_genes, genes_in_object)
perk_genes_filtered <- filter_genes(perk_genes, genes_in_object)
atf6_genes_filtered <- filter_genes(atf6_genes, genes_in_object)
proliferation_genes_filtered <- filter_genes(proliferation_genes, genes_in_object)
mesenchymal_genes_filtered <- filter_genes(mesenchymal_genes, genes_in_object)
stemness_genes_filtered <- filter_genes(stemness_genes, genes_in_object)
hypoxia_genes_filtered <- filter_genes(hypoxia_genes, genes_in_object)

# Print informative messages about the number of genes found
print("Gene lists have been filtered. Number of genes found:")
print(paste("IRE1:", length(ire1_genes_filtered[[1]]), "/", length(ire1_genes)))
print(paste("PERK:", length(perk_genes_filtered[[1]]), "/", length(perk_genes)))
print(paste("ATF6:", length(atf6_genes_filtered[[1]]), "/", length(atf6_genes)))
print(paste("Proliferation:", length(proliferation_genes_filtered[[1]]), "/", length(proliferation_genes)))
print(paste("Mesenchymal:", length(mesenchymal_genes_filtered[[1]]), "/", length(mesenchymal_genes)))
print(paste("Stemness:", length(stemness_genes_filtered[[1]]), "/", length(stemness_genes)))
print(paste("Hypoxia:", length(hypoxia_genes_filtered[[1]]), "/", length(hypoxia_genes)))


# --- CALCULATE MODULE SCORES (using filtered lists) ---

print("Calculating UPR scores...")
seurat_obj <- AddModuleScore(seurat_obj, features = ire1_genes_filtered, name = "IRE1_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = perk_genes_filtered, name = "PERK_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = atf6_genes_filtered, name = "ATF6_Score")

print("Calculating Functional Axis scores...")
seurat_obj <- AddModuleScore(seurat_obj, features = proliferation_genes_filtered, name = "Proliferation_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = mesenchymal_genes_filtered, name = "Mesenchymal_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = stemness_genes_filtered, name = "Stemness_Score")
seurat_obj <- AddModuleScore(seurat_obj, features = hypoxia_genes_filtered, name = "Hypoxia_Score")
print("All scoring is complete.")

# --- SAVE THE UPDATED DATA ---
print("Saving the updated Seurat object: glioblastoma_with_all_scores.rds")
saveRDS(seurat_obj, "glioblastoma_with_all_scores.rds")

print("Process completed successfully.")