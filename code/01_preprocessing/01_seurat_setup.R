# ============================================================
# ADIM 1: SEURAT NESNESI OLUSTURMA
# ============================================================

# --- Gerekli Paketlerin Kurulumu (Akýllý Kontrol) ---
cat("Gerekli R paketleri kontrol ediliyor...\n")

# CRAN Aynasý Ayarý (HPC'de sorunsuz kurulum için)
options(repos = c(CRAN = "https://cran.rstudio.com/"))

if (!requireNamespace("Seurat", quietly = TRUE)) {
  cat("Seurat paketi bulunamadi, kuruluyor...\n")
  install.packages("Seurat")
}
if (!requireNamespace("tidyverse", quietly = TRUE)) {
  cat("Tidyverse paketi bulunamadi, kuruluyor...\n")
  install.packages("tidyverse")
}
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}
cat("Paket kontrolu tamamlandi.\n\n")


# --- Paketleri Yükle ---
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))

cat("Adim 1: Seurat Kurulumu Baslatildi.\n")

# --- 1. Veri Setini Yükle ---
cat("Sayim matrisi okunuyor: glioblastoma_gene_counts.csv\n")
# `check.names=FALSE` parametresi, Seurat'ýn SRR isimlerini deðiþtirmesini engeller.
count_matrix <- read.csv("glioblastoma_gene_counts.csv", row.names = 1, check.names = FALSE)
cat("Matris boyutlari: ", nrow(count_matrix), " gen x ", ncol(count_matrix), " hucre\n")


# --- 2. Seurat Nesnesini Oluþtur ---
# min.cells = 3: Bir genin, en az 3 hücrede saptanmýþ olmasý gerekir.
# min.features = 200: Bir hücrenin, en az 200 farklý gen içermesi gerekir.
cat("Seurat nesnesi olusturuluyor...\n")
glioblastoma_seurat <- CreateSeuratObject(counts = count_matrix, 
                                          project = "GlioblastomaUPR", 
                                          min.cells = 3, 
                                          min.features = 200)

# Oluþturulan nesnenin özetine bakalým
print(glioblastoma_seurat)


# --- 3. Seurat Nesnesini Kaydet ---
# Bu .rds dosyasý, sonraki tüm adýmlarda kullanacaðýmýz ana veri dosyamýz olacak.
cat("Seurat nesnesi kaydediliyor: glioblastoma_seurat.rds\n")
saveRDS(glioblastoma_seurat, file = "glioblastoma_seurat.rds")

cat("\n?? ADIM 1 TAMAMLANDI! ??\n")
cat("Olusturulan Seurat nesnesi 'glioblastoma_seurat.rds' dosyasina kaydedildi.\n")
cat("Sonraki adim: 02_seurat_qc.R ile Kalite Kontrolu.\n")