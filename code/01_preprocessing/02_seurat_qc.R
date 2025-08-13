# ============================================================
# ADIM 2: SEURAT KALITE KONTROLU VE FILTRELEME
# ============================================================

# --- Gerekli Paketleri Yukle ---
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

cat("Adim 2: Seurat Kalite Kontrolu Baslatildi.\n")

# --- 1. Bir onceki adimdaki Seurat nesnesini yukle ---
cat("Seurat nesnesi okunuyor: glioblastoma_seurat.rds\n")
if (!file.exists("glioblastoma_seurat.rds")) {
  stop("HATA: glioblastoma_seurat.rds dosyasi bulunamadi. Lutfen once 01_seurat_setup.R betigini calistirin.")
}
glioblastoma_seurat <- readRDS("glioblastoma_seurat.rds")

# --- 2. Mitokondriyal Gen Yuzdesini Hesapla ---
# Insan mitokondriyal genleri genellikle "MT-" on ekiyle baslar.
cat("Mitokondriyal gen yuzdesi hesaplaniyor...\n")
glioblastoma_seurat[["percent.mt"]] <- PercentageFeatureSet(glioblastoma_seurat, pattern = "^MT-")

# --- 3. Kalite Metriklerini Gorsellestir (Filtrelemeden Once) ---
cat("Kalite metrikleri icin violin grafikleri olusturuluyor...\n")
# Bu PDF dosyasi, filtreleme esiklerinizi belirlemenize yardimci olacaktir.
pdf("QC_metrics_before_filtering.pdf", width = 12, height = 6)
VlnPlot(glioblastoma_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

plot1 <- FeatureScatter(glioblastoma_seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(glioblastoma_seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf("QC_scatter_plots_before_filtering.pdf", width = 12, height = 6)
print(plot1 + plot2)
dev.off()
cat("Gorsellestirmeler 'QC_metrics_before_filtering.pdf' ve 'QC_scatter_plots_before_filtering.pdf' dosyalarina kaydedildi.\n")


# --- 4. Hucreleri Filtrele ---
# NOT: Asagidaki esik degerleri genel baslangic noktalaridir.
# Yukarida olusturulan PDF grafiklerini inceleyerek bu degerleri kendi verinize gore ayarlamaniz ONEMLIDIR!
# nFeature_RNA < 5000: Olasi cift hucreleri (doublet) temizler.
# percent.mt < 10: Stresli veya olmekte olan hucreleri (yuksek mitokondriyal icerik) temizler.
cat("Hucreler filtreleniyor...\n")
cat("Orijinal hucre sayisi: ", ncol(glioblastoma_seurat), "\n")

glioblastoma_filtered <- subset(glioblastoma_seurat, subset = nFeature_RNA > 200 & percent.mt < 40)

cat("Filtrelenmis hucre sayisi: ", ncol(glioblastoma_filtered), "\n")

# --- 5. Filtrelenmis Seurat Nesnesini Kaydet ---
cat("Filtrelenmis Seurat nesnesi kaydediliyor: glioblastoma_filtered.rds\n")
saveRDS(glioblastoma_filtered, file = "glioblastoma_filtered.rds")

cat("\nADIM 2 TAMAMLANDI!\n")
cat("Filtrelenmis Seurat nesnesi 'glioblastoma_filtered.rds' dosyasina kaydedildi.\n")
cat("Sonraki adim: 03_seurat_analysis.R ile Veri Normalizasyonu ve Boyut Indirgeme.\n")