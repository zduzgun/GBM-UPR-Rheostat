# ==============================================================================
# 07_01_neftel_data_loading_and_setup.R (v4 - Alternative Metadata Approach)
#
# Purpose: To load the Neftel et al. (2019) 10x GBM data and create a Seurat object
#          This version handles cases where proper metadata might not be available
# ==============================================================================

# --- 1. Load Libraries ---
cat("Step 1: Loading libraries...\n")
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readxl)
  library(data.table)
  library(stringr)
})

# --- 2. Define File Paths ---
cat("Step 2: Defining file paths...\n")
data_file <- "GSM3828673_10X_GBM_IDHwt_processed_TPM.tsv"
metadata_file <- "GSE131928_single_cells_tumor_name_and_adult_or_peidatric.xlsx"
output_rds_file <- "07_01_neftel_initial_seurat_object.rds"

# --- 3. Load Data ---
cat("Step 3: Loading gene expression matrix (TPM)...\n")
neftel_tpm <- fread(data_file, data.table = FALSE)

# Set gene names as row names and remove the first column
rownames(neftel_tpm) <- neftel_tpm[, 1]
neftel_tpm <- neftel_tpm[, -1]

cat("Data matrix dimensions:", dim(neftel_tpm)[1], "genes x", dim(neftel_tpm)[2], "cells\n")

# --- 4. Try to Load Metadata (with fallback options) ---
cat("Step 4: Attempting to load metadata...\n")

# Function to try different approaches to read the metadata
try_read_metadata <- function(file_path) {
  approaches <- list(
    # Approach 1: Try different sheet names
    function() {
      tryCatch({
        sheets <- excel_sheets(file_path)
        cat("Available Excel sheets:", paste(sheets, collapse = ", "), "\n")
        
        # Look for sheets that might contain actual data
        data_sheet <- NULL
        for (sheet in sheets) {
          if (grepl("data|sample|cell|metadata", sheet, ignore.case = TRUE)) {
            data_sheet <- sheet
            break
          }
        }
        
        if (!is.null(data_sheet)) {
          cat("Trying to read sheet:", data_sheet, "\n")
          return(read_excel(file_path, sheet = data_sheet))
        }
        
        # If no obvious sheet, try the last one (often contains data)
        if (length(sheets) > 1) {
          last_sheet <- sheets[length(sheets)]
          cat("Trying last sheet:", last_sheet, "\n")
          return(read_excel(file_path, sheet = last_sheet))
        }
        
        return(NULL)
      }, error = function(e) NULL)
    },
    
    # Approach 2: Skip many rows to find actual data
    function() {
      tryCatch({
        for (skip_rows in c(10, 20, 50, 100, 200)) {
          cat("Trying to skip", skip_rows, "rows...\n")
          temp_data <- read_excel(file_path, skip = skip_rows)
          
          # Check if this looks like real data
          if (ncol(temp_data) > 5 && 
              any(grepl("cell|sample|barcode", colnames(temp_data), ignore.case = TRUE))) {
            cat("Found potential data starting at row", skip_rows + 1, "\n")
            return(temp_data)
          }
        }
        return(NULL)
      }, error = function(e) NULL)
    },
    
    # Approach 3: Look for CSV version
    function() {
      csv_file <- gsub("\\.xlsx?$", ".csv", file_path)
      if (file.exists(csv_file)) {
        cat("Found CSV version, trying to read:", csv_file, "\n")
        return(read.csv(csv_file))
      }
      return(NULL)
    }
  )
  
  for (approach in approaches) {
    result <- approach()
    if (!is.null(result) && nrow(result) > 0) {
      return(result)
    }
  }
  
  return(NULL)
}

# Try to read metadata
meta <- try_read_metadata(metadata_file)

if (is.null(meta)) {
  cat("WARNING: Could not read proper metadata from the Excel file.\n")
  cat("Creating minimal metadata from cell names...\n")
  
  # Create basic metadata from cell barcodes
  cell_names <- colnames(neftel_tpm)
  meta <- data.frame(
    Cell.name = cell_names,
    Sample = "Neftel_GBM",
    Technology = "10x",
    row.names = cell_names,
    stringsAsFactors = FALSE
  )
  
  # Try to extract additional info from cell names if possible
  if (any(grepl("_", cell_names))) {
    # If cell names have underscores, try to extract sample info
    sample_info <- str_extract(cell_names, "^[^_]+")
    if (length(unique(sample_info)) > 1) {
      meta$Sample <- sample_info
    }
  }
  
} else {
  cat("Successfully loaded metadata with", nrow(meta), "rows and", ncol(meta), "columns\n")
  cat("Column names:", paste(colnames(meta), collapse = ", "), "\n")
  
  # Clean column names
  colnames(meta) <- make.names(colnames(meta), unique = TRUE)
  
  # Try to identify the cell name column
  cell_col <- NULL
  possible_cell_cols <- c("Cell.name", "cell.name", "Cell", "cell", "Barcode", "barcode", 
                         "Sample.name", "sample.name")
  
  for (col in possible_cell_cols) {
    if (col %in% colnames(meta)) {
      cell_col <- col
      break
    }
  }
  
  if (is.null(cell_col)) {
    # Look for columns that might contain cell identifiers
    for (i in 1:ncol(meta)) {
      col_values <- as.character(meta[, i])
      # Check if values look like cell barcodes or names
      if (sum(col_values %in% colnames(neftel_tpm)) > 100) {
        cell_col <- colnames(meta)[i]
        break
      }
    }
  }
  
  if (!is.null(cell_col)) {
    cat("Using column '", cell_col, "' as cell identifier\n")
    
    # Filter for cells that exist in expression data
    valid_cells <- meta[[cell_col]] %in% colnames(neftel_tpm)
    meta <- meta[valid_cells, ]
    rownames(meta) <- meta[[cell_col]]
    
    cat("Found", nrow(meta), "cells with matching expression data\n")
  } else {
    cat("WARNING: Could not identify cell name column in metadata\n")
    cat("Creating basic metadata instead...\n")
    
    # Fallback to basic metadata
    cell_names <- colnames(neftel_tpm)
    meta <- data.frame(
      Cell.name = cell_names,
      Sample = "Neftel_GBM",
      Technology = "10x",
      row.names = cell_names,
      stringsAsFactors = FALSE
    )
  }
}

# Ensure we have matching cells
common_cells <- intersect(colnames(neftel_tpm), rownames(meta))
cat("Common cells between expression data and metadata:", length(common_cells), "\n")

if (length(common_cells) == 0) {
  cat("WARNING: No matching cells found. Using all expression data cells.\n")
  cell_names <- colnames(neftel_tpm)
  meta <- data.frame(
    Cell.name = cell_names,
    Sample = "Neftel_GBM",
    Technology = "10x",
    row.names = cell_names,
    stringsAsFactors = FALSE
  )
  common_cells <- cell_names
}

# Subset data to common cells
neftel_tpm_filtered <- neftel_tpm[, common_cells]
meta_filtered <- meta[common_cells, , drop = FALSE]

cat("Final dimensions - Expression:", dim(neftel_tpm_filtered), ", Metadata:", dim(meta_filtered), "\n")

# --- 5. Create Seurat Object ---
cat("Step 5: Creating Seurat object...\n")
neftel_seurat <- CreateSeuratObject(counts = neftel_tpm_filtered,
                                    project = "Neftel_Validation",
                                    meta.data = meta_filtered)

cat("Seurat object created with", ncol(neftel_seurat), "cells and", nrow(neftel_seurat), "genes.\n")

# --- 6. Basic Preprocessing and QC ---
cat("Step 6: Applying Quality Control (QC) and preprocessing steps...\n")
DefaultAssay(neftel_seurat) <- "RNA"

# Add mitochondrial gene percentage
neftel_seurat[["percent.mt"]] <- PercentageFeatureSet(neftel_seurat, pattern = "^MT-")

# Add ribosomal gene percentage
neftel_seurat[["percent.ribo"]] <- PercentageFeatureSet(neftel_seurat, pattern = "^RP[SL]")

# Calculate basic QC metrics
neftel_seurat[["nFeature_RNA"]] <- nrow(neftel_seurat)
neftel_seurat[["nCount_RNA"]] <- colSums(neftel_tpm_filtered)

# Normalize data
cat("   - Normalizing data...\n")
neftel_seurat <- NormalizeData(neftel_seurat, normalization.method = "LogNormalize", scale.factor = 1e6)

# Find variable features
cat("   - Finding variable features...\n")
neftel_seurat <- FindVariableFeatures(neftel_seurat, selection.method = "vst", nfeatures = 2000)

cat("   - Scaling data (this may take a while)...\n")
all.genes <- rownames(neftel_seurat)
neftel_seurat <- ScaleData(neftel_seurat, features = all.genes)

cat("   - Running PCA...\n")
neftel_seurat <- RunPCA(neftel_seurat, features = VariableFeatures(object = neftel_seurat))

# --- 7. Save the Result ---
cat("Step 7: Saving the processed Seurat object...\n")
saveRDS(neftel_seurat, file = output_rds_file)

cat("Analysis complete. Output saved to '", output_rds_file, "'.\n")
cat("Summary:\n")
cat("  - Total cells:", ncol(neftel_seurat), "\n")
cat("  - Total genes:", nrow(neftel_seurat), "\n")
cat("  - Variable features:", length(VariableFeatures(neftel_seurat)), "\n")
cat("  - Available metadata columns:", paste(colnames(meta_filtered), collapse = ", "), "\n")

# Print some basic statistics
cat("\nBasic QC statistics:\n")
cat("  - Median genes per cell:", median(neftel_seurat$nFeature_RNA), "\n")
cat("  - Median counts per cell:", median(neftel_seurat$nCount_RNA), "\n")
cat("  - Median mitochondrial %:", median(neftel_seurat$percent.mt), "\n")