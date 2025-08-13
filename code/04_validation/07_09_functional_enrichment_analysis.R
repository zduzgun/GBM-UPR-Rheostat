# ==============================================================================
# 07_09_functional_enrichment_analysis.R (v2 - Dedicated Output Folder)
#
# Purpose: Perform functional enrichment analysis (GO and KEGG) and save all
#          outputs into a single, dedicated folder for better organization.
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db) # Human-specific gene annotation database
  library(dplyr)
  library(ggplot2)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
# --- GÜNCELLEME BURADA: Çýktý klasörünü ve yollarýný deðiþtir ---
output_directory <- "07_09_functional_enrichment_results"
input_gene_list_csv <- "07_08_pseudotime_gene_analysis_results/tables/07_08_pseudotime_dependent_genes.csv"
output_go_results_csv <- file.path(output_directory, "07_09_GO_enrichment_results.csv")
output_kegg_results_csv <- file.path(output_directory, "07_09_KEGG_enrichment_results.csv")
output_plot_file <- file.path(output_directory, "07_09_enrichment_plots.pdf")

# --- GÜNCELLEME BURADA: Yeni çýktý klasörünü oluþtur ---
# Çýktý klasörü mevcut deðilse, betiðin baþýnda oluþtur.
cat("Step 2.1: Ensuring output directory exists...\n")
dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)


# --- 3. Load and Prepare the Gene List ---
cat("Step 3: Loading and preparing the list of significant genes...\n")
all_sig_genes <- read.csv(input_gene_list_csv)
cat("  - Loaded", nrow(all_sig_genes), "significant genes from the previous step.\n")

gene_symbols <- all_sig_genes$gene_short_name
entrez_ids <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
cat("  - Successfully converted", nrow(entrez_ids), "gene symbols to ENTREZ IDs.\n")


# --- 4. Perform GO (Gene Ontology) Enrichment Analysis ---
cat("Step 4: Performing GO Enrichment Analysis...\n")
go_results <- enrichGO(gene          = entrez_ids$ENTREZID,
                         OrgDb         = org.Hs.eg.db,
                         ont           = "ALL",
                         qvalueCutoff  = 0.05,
                         readable      = TRUE)

if (!is.null(go_results)) {
  write.csv(as.data.frame(go_results), file = output_go_results_csv, row.names = FALSE)
  cat("  - GO analysis complete. Found", nrow(go_results), "enriched terms.\n")
} else {
  cat("  - No significant GO terms found.\n")
}


# --- 5. Perform KEGG Pathway Enrichment Analysis ---
cat("Step 5: Performing KEGG Pathway Enrichment Analysis...\n")
kegg_results <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                           organism     = 'hsa',
                           qvalueCutoff  = 0.05)

if (!is.null(kegg_results)) {
  kegg_results <- setReadable(kegg_results, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  write.csv(as.data.frame(kegg_results), file = output_kegg_results_csv, row.names = FALSE)
  cat("  - KEGG analysis complete. Found", nrow(kegg_results), "enriched pathways.\n")
} else {
  cat("  - No significant KEGG pathways found.\n")
}


# --- 6. Visualize Enrichment Results ---
cat("Step 6: Generating visualizations...\n")
pdf(output_plot_file, width = 14, height = 10)

if (!is.null(go_results) && nrow(go_results) > 0) {
  p_go <- dotplot(go_results, showCategory = 20, split="ONTOLOGY") +
          facet_grid(ONTOLOGY~., scale="free") +
          ggtitle("Top Enriched Gene Ontology Terms")
  print(p_go)
}

if (!is.null(kegg_results) && nrow(kegg_results) > 0) {
  p_kegg <- dotplot(kegg_results, showCategory = 20) +
            ggtitle("Top Enriched KEGG Pathways")
  print(p_kegg)
}

dev.off()
cat("Analysis complete. Outputs are saved in:", output_directory, "\n")