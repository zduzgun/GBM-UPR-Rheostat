# =================================================================
# ADIM 4: UPR REOSTAT HIPOTEZININ ANALIZI
# =================================================================

# --- Gerekli Paketleri Yukle ---
suppressMessages(library(Seurat))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

cat("Adim 4: UPR Reostat Analizi Baslatildi.\n")

# --- 1. Islenmis Seurat nesnesini yukle ---
cat("Islenmis Seurat nesnesi okunuyor: glioblastoma_processed.rds\n")
if (!file.exists("glioblastoma_processed.rds")) {
  stop("HATA: glioblastoma_processed.rds dosyasi bulunamadi. Lutfen once 03_seurat_analysis.R betigini calistirin.")
}
glio.seu <- readRDS("glioblastoma_processed.rds")

# --- 2. UPR Kollari icin Gen Listeleri Olustur ---
# NOT: Bu gen listeleri, projenizin bilimsel temelini olusturur.
# Yayina donuk bir calisma icin bu listelerin GSEA/MSigDB gibi veritabanlarindan
# veya kapsamli literatur taramasindan dikkatlice kurate edilmesi gerekir.
# Buradakiler, yola baslamak icin ornek listelerdir.
cat("UPR yolagi gen listeleri olusturuluyor...\n")

# IRE1? yolu (XBP1s hedefleri)
ire1_genes <- list(c("XBP1", "ERN1", "DNAJB9", "EDEM1", "HSPA5", "PDIA4"))

# PERK yolu (ATF4 hedefleri)
perk_genes <- list(c("EIF2AK3", "ATF4", "DDIT3", "TRIB3", "ASNS", "PSAT1"))

# ATF6 yolu (ATF6 hedefleri)
atf6_genes <- list(c("ATF6", "CALR", "MANF", "HERPUD1", "HSP90B1"))

# --- 3. Her Hucre Icin Yolagi Aktivite Skorlarini Hesapla ---
# AddModuleScore fonksiyonu, belirtilen genlerin ortalama ifadesini,
# rastgele secilmis kontrol genlerine gore duzelterek bir skor hesaplar.
cat("Her bir UPR kolu icin aktivite skorlari hesaplaniyor...\n")

glio.seu <- AddModuleScore(object = glio.seu, features = ire1_genes, name = 'IRE1_Score')
glio.seu <- AddModuleScore(object = glio.seu, features = perk_genes, name = 'PERK_Score')
glio.seu <- AddModuleScore(object = glio.seu, features = atf6_genes, name = 'ATF6_Score')

# --- 4. Skorlari UMAP uzerinde Gorsellestir ---
# Bu gorsel, reostat hipotezinin en dogrudan kaniti olacaktir.
cat("UPR skorlari UMAP uzerinde gorsellestiriliyor...\n")

p1 <- FeaturePlot(glio.seu, features = "IRE1_Score1", pt.size = 0.5, order = TRUE) + ggtitle("IRE1 Yolagi Aktivitesi")
p2 <- FeaturePlot(glio.seu, features = "PERK_Score1", pt.size = 0.5, order = TRUE) + ggtitle("PERK Yolagi Aktivitesi")
p3 <- FeaturePlot(glio.seu, features = "ATF6_Score1", pt.size = 0.5, order = TRUE) + ggtitle("ATF6 Yolagi Aktivitesi")

# Grafikleri tek bir dosyaya kaydet
pdf("UPR_Reostat_UMAP_Visualization.pdf", width = 18, height = 6)
print(p1 | p2 | p3)
dev.off()
cat("Gorsellestirmeler 'UPR_Reostat_UMAP_Visualization.pdf' dosyasina kaydedildi.\n")

# --- 5. Skorlarin Dagilimini Kumelere Gore Incele (Violin Plots) ---
cat("Skorlarin kumelere gore dagilimi icin violin grafikleri olusturuluyor...\n")
pdf("UPR_Reostat_VlnPlots.pdf", width = 12, height = 5)
VlnPlot(glio.seu, features = c("IRE1_Score1", "PERK_Score1", "ATF6_Score1"), pt.size = 0, ncol = 3)
dev.off()
cat("Violin grafikleri 'UPR_Reostat_VlnPlots.pdf' dosyasina kaydedildi.\n")

cat("\n?? ADIM 4 TAMAMLANDI! ??\n")
cat("UPR Reostat analizi tamamlandi.\n")
cat("Sonuclari 'UPR_Reostat_UMAP_Visualization.pdf' ve 'UPR_Reostat_VlnPlots.pdf' dosyalarinda bulabilirsiniz.\n")