# =====================================================================
# ADIM 5: HÜCRE KÜMELERÝNÝN KARAKTERÝZASYONU VE MARKER ANALÝZÝ
# =====================================================================

# --- Gerekli Paketleri Yukle ---
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

cat("Adim 5: Kume Karakterizasyonu Baslatildi.\n")

# --- 1. Bir onceki adimdaki islenmis Seurat nesnesini yukle ---
cat("Islenmis Seurat nesnesi okunuyor: glioblastoma_processed.rds\n")
if (!file.exists("glioblastoma_processed.rds")) {
  stop("HATA: glioblastoma_processed.rds dosyasi bulunamadi. Lutfen once 03_seurat_analysis.R betigini calistirin.")
}
glio.seu <- readRDS("glioblastoma_processed.rds")

# --- 2. Her kume icin marker genleri bul ---
# Bu islem, cekirdek sayisina ve veri buyuklugune gore biraz zaman alabilir.
cat("Her bir kume icin marker genler bulunuyor (FindAllMarkers)...\n")
# Sadece pozitif marker'lari buluyoruz ve en az %25 hucrede ifade edilenleri aliyoruz.
all_markers <- FindAllMarkers(glio.seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
cat("Marker genler bulundu.\n")

# --- 3. Marker genleri CSV dosyasina kaydet ---
# Bu dosya, hucre tipi tanimlamasi icin en onemli ciktidir.
cat("Marker genler 'all_cluster_markers.csv' dosyasina kaydediliyor...\n")
write.csv(all_markers, file = "all_cluster_markers.csv", row.names = FALSE)

# --- 4. Gorsellestirme icin en iyi marker'lari sec ---
# Her kumeden, ifade farki en yuksek olan ilk 10 marker geni aliyoruz.
top10_markers <- all_markers %>%
    group_by(cluster) %>%
    slice_max(n = 10, order_by = avg_log2FC)

# --- 5. Heatmap (Isi Haritasi) olustur ---
cat("En iyi marker genler icin heatmap olusturuluyor...\n")
pdf("marker_heatmap.pdf", width = 16, height = 12)
DoHeatmap(glio.seu, features = top10_markers$gene) + NoLegend()
dev.off()
cat("Heatmap 'marker_heatmap.pdf' dosyasina kaydedildi.\n")

# --- 6. Dot Plot olustur ---
# Dot plot, gen ifadesinin kumelerdeki yayginligini (nokta boyutu) ve ortalama ifadesini (renk) gosterir.
cat("En iyi marker genler icin dot plot olusturuluyor...\n")
pdf("marker_dotplot.pdf", width = 16, height = 8)
# unique() kullanarak, bir genin birden fazla kume icin top10 olmasi durumunda tekrari onluyoruz.
DotPlot(glio.seu, features = unique(top10_markers$gene)) + RotatedAxis()
dev.off()
cat("Dot plot 'marker_dotplot.pdf' dosyasina kaydedildi.\n")


cat("\n?? ADIM 5 TAMAMLANDI! ??\n")
cat("Marker genler ve gorseller olusturuldu.\n")
cat("Simdi 'all_cluster_markers.csv' dosyasini inceleyerek kumelere biyolojik kimlik atayabilirsiniz.\n")