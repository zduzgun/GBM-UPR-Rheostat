# =====================================================================
# ADIM 3: VERI NORMALIZASYONU, BOYUT INDIRGEME VE KUMELEME
# =====================================================================

# --- Gerekli Paketleri Yukle ---
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

cat("Adim 3: Veri Analizi Baslatildi.\n")

# --- 1. Filtrelenmis Seurat nesnesini yukle ---
cat("Filtrelenmis Seurat nesnesi okunuyor: glioblastoma_filtered.rds\n")
if (!file.exists("glioblastoma_filtered.rds")) {
  stop("HATA: glioblastoma_filtered.rds dosyasi bulunamadi. Lutfen once 02_seurat_qc.R betigini calistirin.")
}
glio.seu <- readRDS("glioblastoma_filtered.rds")

# --- 2. Veriyi Normalize Et ---
# Bu adim, hucreler arasi sekanslama derinligi farkliliklarini hesaba katar.
cat("Veri normalizasyonu yapiliyor (LogNormalize)...\n")
glio.seu <- NormalizeData(glio.seu)

# --- 3. Yuksek Varyansli Genleri Bul ---
# Hucreler arasi varyasyonu en cok aciklayan genleri belirler. Bu genler analizin temelini olusturur.
cat("Yuksek varyansli genler bulunuyor...\n")
glio.seu <- FindVariableFeatures(glio.seu, selection.method = "vst", nfeatures = 2000)

# --- 4. Veriyi Olceklendir (Scaling) ---
# Bu adim, yuksek ifade edilen genlerin analizi domine etmesini engeller. Sadece degisken genler olceklendirilir.
cat("Veri olceklendiriliyor...\n")
all.genes <- rownames(glio.seu)
glio.seu <- ScaleData(glio.seu, features = all.genes)

# --- 5. Lineer Boyut Indirgeme (PCA) ---
# Binlerce gen boyutunu, verinin varyasyonunu en cok aciklayan daha az sayida temel bilesene (PC) indirger.
cat("PCA calistiriliyor...\n")
glio.seu <- RunPCA(glio.seu, features = VariableFeatures(object = glio.seu))

# --- 6. Kac adet PC kullanilacagina karar ver (Elbow Plot) ---
# Genellikle grafigin "dirsek" yaptigi nokta secilir. Bu, sinyalin gurultuden ayrildigi yerdir.
cat("Elbow plot olusturuluyor...\n")
pdf("elbow_plot.pdf")
ElbowPlot(glio.seu)
dev.off()
cat("Gorsel 'elbow_plot.pdf' dosyasina kaydedildi.\n")

# --- 7. Hucreleri Kumele ---
# PCA sonuclarini kullanarak hucreler arasi mesafeleri hesaplar ve hucreleri kumelere ayirir.
# NOT: `dims` parametresini elbow plot'a bakarak ayarlayabilirsiniz. Genellikle 15-30 arasi bir deger kullanilir.
cat("Hucre kumeleri bulunuyor...\n")
glio.seu <- FindNeighbors(glio.seu, dims = 1:20)
glio.seu <- FindClusters(glio.seu, resolution = 0.5)

# --- 8. Non-lineer Boyut Indirgeme (UMAP) ---
# Kumeleri 2 boyutlu bir harita uzerinde gorsellestirmek icin kullanilir.
cat("UMAP calistiriliyor...\n")
glio.seu <- RunUMAP(glio.seu, dims = 1:20)

# --- 9. UMAP Gorselini Kaydet ---
cat("UMAP gorseli olusturuluyor...\n")
pdf("UMAP_clusters.pdf", width = 8, height = 7)
DimPlot(glio.seu, reduction = "umap", label = TRUE)
dev.off()
cat("Gorsel 'UMAP_clusters.pdf' dosyasina kaydedildi.\n")

# --- 10. Islenmis Seurat Nesnesini Kaydet ---
cat("Islenmis Seurat nesnesi kaydediliyor: glioblastoma_processed.rds\n")
saveRDS(glio.seu, file = "glioblastoma_processed.rds")

cat("\n?? ADIM 3 TAMAMLANDI! ??\n")
cat("Islenmis Seurat nesnesi 'glioblastoma_processed.rds' dosyasina kaydedildi.\n")
cat("Sonraki adim: 04_upr_reostat_analysis.R ile hipotez testi.\n")