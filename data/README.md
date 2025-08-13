# Data Files

All data files for this project are available on Zenodo due to size constraints.

## Download Instructions

### Option 1: Download Complete Dataset
```bash
# Download from Zenodo (DOI will be provided after publication)
wget https://zenodo.org/record/XXXXX/files/GBM-UPR-Rheostat-Data.tar.gz
tar -xzf GBM-UPR-Rheostat-Data.tar.gz
```

### Option 2: Download Specific Files
Individual files can be downloaded from the Zenodo repository web interface.

## Data Structure

```
data/
+¦¦ processed/      # 17 GB - Seurat and analysis objects
+¦¦ raw/           # 1.8 GB - Original data from GEO
+¦¦ referans/      # 4.3 GB - Reference genomes and indices
L¦¦ resources/     # 392 MB - Gene lists and databases
```

## File Descriptions

### Key Processed Files
- `glioblastoma_seurat.rds` (106 MB) - Initial Seurat object
- `glioblastoma_with_final_scores.rds` (468 MB) - Final analysis with UPR scores
- `07_03_neftel_with_scores.rds` (3.0 GB) - Validation dataset with scores

### Raw Data Sources
- GSE57872 - Patel et al. 2014 dataset
- GSE131928 - Neftel et al. 2019 dataset

## Citation
If using this data, please cite both the original studies and our analysis.
