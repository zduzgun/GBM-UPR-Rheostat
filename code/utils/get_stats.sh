#!/bin/bash

#-----------------------------------------------------------------------
# HEDEFE Y�NEL�K METODOLOJ� B�LG� TOPLAMA BET��� v2.0
# Bu betik, sadece ilgili betik dosyalar�n� tarar ve anahtar yaz�l�m
# versiyonlar�n� listeleyerek k���k ve �zet bir log dosyas� olu�turur.
#-----------------------------------------------------------------------

LOG_FILE="statistic_summary.log"

# Log dosyas�n� ba�lat ve tarih damgas� ekle
echo "Methodological Details Summary - Generated on $(date)" > "$LOG_FILE"
echo "=======================================================" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# --- 1. ANAHTAR ANAL�Z PARAMETRELER� ---
echo "--- 1. Key Analysis Parameters (from .R scripts) ---" >> "$LOG_FILE"
echo "[Patel] PCA Dims & UMAP (03_seurat_analysis.R):" >> "$LOG_FILE"
grep -E "FindNeighbors|RunUMAP" 03_seurat_analysis.R >> "$LOG_FILE" 2>&1
echo "" >> "$LOG_FILE"
echo "[Patel] Clustering Resolution (03_seurat_analysis.R):" >> "$LOG_FILE"
grep "FindClusters" 03_seurat_analysis.R >> "$LOG_FILE" 2>&1
echo "" >> "$LOG_FILE"
echo "[Neftel] Initial & Malignant PCA Dims (relevant .R files):" >> "$LOG_FILE"
find ./07_Validation_Neftel -name "*.R" -print0 | xargs -0 grep -H "RunUMAP" >> "$LOG_FILE" 2>&1
echo "" >> "$LOG_FILE"
echo "[Seurat] Normalization Method Check (All .R files):" >> "$LOG_FILE"
find . -name "*.R" -print0 | xargs -0 grep -H -E "NormalizeData|SCTransform" >> "$LOG_FILE" 2>&1
echo "" >> "$LOG_FILE"
echo "[Reproducibility] Random Seed Check (All .R files):" >> "$LOG_FILE"
find . -name "*.R" -print0 | xargs -0 grep -H "set.seed" >> "$LOG_FILE" 2>&1
echo "-------------------------------------------------------" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# --- 2. KALLISTO & BATCH CORRECTION KONTROL� ---
echo "--- 2. Kallisto & Batch Correction Parameters ---" >> "$LOG_FILE"
echo "[Kallisto] Bootstrap parameter (-b) in run_kallisto.py:" >> "$LOG_FILE"
grep "kallisto quant" run_kallisto.py >> "$LOG_FILE" 2>&1
echo "" >> "$LOG_FILE"
echo "[QC] Batch Correction Check (in .R files only):" >> "$LOG_FILE"
find . -name "*.R" -print0 | xargs -0 grep -H -i -E "harmony|integrate|cca|liger" >> "$LOG_FILE" 2>&1
echo "-------------------------------------------------------" >> "$LOG_FILE"
echo "" >> "$LOG_FILE"

# --- 3. ANAHTAR YAZILIM VERS�YONLARI ---
echo "--- 3. Key Software Versions (Filtered from Conda) ---" >> "$LOG_FILE"
KEY_PACKAGES="r-base|seurat|monocle|cellchat|pyscenic|kallisto|sra-tools|fastqc|multiqc|python|tximport"
ENVS=("bioinfo_env" "monocle_env" "scenic_env" "seurat_env")
for ENV_NAME in "${ENVS[@]}"; do
    echo "[Conda] Key packages in: $ENV_NAME" >> "$LOG_FILE"
    # 'conda list' komutu bazen hata verebilir, bu y�zden 2>/dev/null ekliyoruz
    conda list -n "$ENV_NAME" 2>/dev/null | grep -i -E "$KEY_PACKAGES" >> "$LOG_FILE"
    echo "" >> "$LOG_FILE"
done
echo "-------------------------------------------------------" >> "$LOG_FILE"

echo "Targeted script finished. Summary of results saved to $LOG_FILE"