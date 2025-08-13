# The Dynamic UPR Rheostat Orchestrates Single-Cell Plasticity in Glioblastoma

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/PENDING.svg)](https://doi.org/PENDING)

## Overview

This repository contains the complete computational workflow for single-cell RNA-seq analysis investigating the Unfolded Protein Response (UPR) rheostat in glioblastoma. Our analysis reveals arm-resolved UPR control mechanisms that orchestrate cellular plasticity in GBM.

## Key Findings

- The UPR operates as a graded, arm-resolved rheostat rather than a binary switch
- IRE1/XBP1 and ATF6 arms are tightly coordinated while PERK exhibits semi-independent behavior
- Rheostatic tuning correlates with hypoxia and metabolic programs
- Cross-dataset validation confirms robustness across technologies (SMART-seq and 10x)

## Data

- **Discovery Cohort**: Patel et al. 2014 (GSE57872) - 871 malignant cells
- **Validation Cohort**: Neftel et al. 2019 (GSE131928) - 11,877 malignant cells

## Repository Structure

```
GBM-UPR-Rheostat/
+¦¦ code/               # Analysis scripts
+¦¦ data/               # Processed data and resources
+¦¦ results/            # Figures and tables
+¦¦ environment/        # Package dependencies
L¦¦ docs/               # Documentation
```

## Requirements

- R >= 4.3.0
- Seurat >= 5.3.0
- Monocle3 >= 1.4.0
- Additional packages listed in `environment/requirements.txt`

## Usage

Scripts are numbered in order of execution:

1. Preprocessing: `code/01_preprocessing/`
2. Main Analysis: `code/02_main_analysis/`
3. Integrated Analysis: `code/03_integrated_analysis/`
4. Validation: `code/04_validation/`
5. Figure Generation: `code/05_figures/`

## Citation

If you use this code or data, please cite:

```
Duzgun Z. (2025). The Dynamic UPR Rheostat Orchestrates Single-Cell Plasticity in Glioblastoma. 
[Journal Name - Pending]. DOI: [Pending]
```

## Author

**Zekeriya Duzgun**  
Department of Medical Biology, Faculty of Medicine  
Giresun University, Giresun, Türkiye  
Email: zekeriya.duzgun@giresun.edu.tr, zduzgun@gmail.com

## Acknowledgments

Computational resources were provided by TUBITAK ULAKBIM, High Performance and Grid Computing Center (TRUBA).

## License

This project is licensed under the MIT License - see the LICENSE file for details.
