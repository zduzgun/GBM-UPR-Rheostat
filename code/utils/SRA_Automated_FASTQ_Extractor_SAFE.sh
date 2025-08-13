#!/bin/bash

# ============================================================
# GNU PARALLEL İLE GÜÇLENDİRİLMİŞ SRA->FASTQ DÖNÜŞTÜRME
# İYİLEŞTİRİLMİŞ VERSİYON
# ============================================================

# --- YAPILANDIRMA DEĞİŞKENLERİ ---
OUT_DIR="fastq_files"
ACC_LIST="resmi_liste_875.txt"
LOG_FILE="hatali_dosyalar.log"
PARALLEL_LOG="parallel_job.log"

# Paralel işlem sayısı (CPU çekirdek sayısı kadar)
MAX_PARALLEL=$(nproc)
MAX_RETRIES=3
MIN_FILE_SIZE=1024  # 1KB minimum dosya boyutu

# --- YARDIMCI FONKSİYONLAR ---

# Temizlik fonksiyonu (sinyal yakalama için)
cleanup() {
    echo -e "\n\n🛑 İşlem durduruldu. Temizlik yapılıyor..."
    
    # Aktif parallel işlemlerini durdur
    pkill -f "parallel.*process_srr" 2>/dev/null
    
    # Yarım kalmış dosyaları temizle
    find "$OUT_DIR" -name "*.tmp" -delete 2>/dev/null
    find "$OUT_DIR" -size -${MIN_FILE_SIZE}c -name "*.fastq.gz" -delete 2>/dev/null
    
    echo "✅ Temizlik tamamlandı. Resume ile devam edebilirsiniz."
    exit 130
}

# Sinyal yakalama
trap cleanup SIGINT SIGTERM

# Resume kontrolü (geliştirilmiş)
check_resume_status() {
    if [ -f "$PARALLEL_LOG" ]; then
        local completed=$(awk 'NR > 1 && $7 == 0' "$PARALLEL_LOG" | wc -l)
        local failed=$(awk 'NR > 1 && $7 != 0' "$PARALLEL_LOG" | wc -l)
        local total=$(wc -l < "$ACC_LIST")
        local remaining=$((total - completed - failed))
        
        if [ "$remaining" -gt 0 ]; then
            echo "🔄 Önceki çalışma tespit edildi:"
            echo "   ✅ Tamamlanan: $completed"
            echo "   ❌ Başarısız: $failed"
            echo "   ⏳ Kalan: $remaining"
            echo "   📊 Resume ile devam edilecek..."
            return 0
        fi
    fi
    return 1
}

# Çıktı dosyası doğrulama (geliştirilmiş)
verify_output() {
    local srr=$1
    local verified=0
    
    # Tüm olası çıktı dosyalarını kontrol et
    for suffix in ".fastq.gz" "_1.fastq.gz" "_2.fastq.gz"; do
        local file="${OUT_DIR}/${srr}${suffix}"
        if [ -f "$file" ]; then
            local size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0)
            if [ "$size" -gt "$MIN_FILE_SIZE" ]; then
                # Gzip formatını kontrol et
                if gzip -t "$file" 2>/dev/null; then
                    verified=$((verified + 1))
                else
                    echo "⚠️ $file bozuk gzip, siliniyor..."
                    rm -f "$file"
                fi
            else
                echo "⚠️ $file çok küçük ($size bytes), siliniyor..."
                rm -f "$file"
            fi
        fi
    done
    
    return $((verified > 0 ? 0 : 1))
}

# Ana işleme fonksiyonu (geliştirilmiş)
process_srr() {
    local srr=$1
    local job_num=$2
    local total_jobs=$3
    
    # Resume kontrolü - dosya zaten varsa atla
    if [ -f "${OUT_DIR}/${srr}_1.fastq.gz" ] || [ -f "${OUT_DIR}/${srr}.fastq.gz" ]; then
        echo "($job_num/$total_jobs) ⏭️ $srr zaten mevcut, atlanıyor."
        return 0
    fi
    
    echo "($job_num/$total_jobs) 🔄 $srr işleniyor..."
    
    # fastq-dump ile timeout koruması
    if timeout 1800 fastq-dump --gzip --split-files --outdir "$OUT_DIR" --origfmt "$srr" 2>/dev/null; then
        # Çıktıyı doğrula
        if verify_output "$srr"; then
            echo "($job_num/$total_jobs) ✅ $srr başarıyla dönüştürüldü."
            return 0
        else
            echo "($job_num/$total_jobs) ❌ $srr doğrulama başarısız."
            return 1
        fi
    else
        echo "($job_num/$total_jobs) ❌ $srr fastq-dump hatası veya timeout."
        # Başarsız durumda yarım kalmış dosyaları temizle
        rm -f "${OUT_DIR}/${srr}"*.fastq.gz "${OUT_DIR}/${srr}"*.tmp 2>/dev/null
        return 1
    fi
}

# Fonksiyonları GNU Parallel için export et
export -f process_srr verify_output
export OUT_DIR MIN_FILE_SIZE

# ============================================================
# ANA PROGRAM
# ============================================================

# Başlangıç kontrolleri
if [ ! -f "$ACC_LIST" ]; then
    echo "❌ Hata: $ACC_LIST dosyası bulunamadı!"
    exit 1
fi

if ! command -v fastq-dump >/dev/null 2>&1; then
    echo "❌ Hata: fastq-dump komutu bulunamadı! SRA Toolkit kurulu mu?"
    exit 1
fi

if ! command -v parallel >/dev/null 2>&1; then
    echo "❌ Hata: GNU parallel komutu bulunamadı!"
    echo "Ubuntu/Debian: sudo apt install parallel"
    echo "CentOS/RHEL: sudo yum install parallel"
    echo "macOS: brew install parallel"
    exit 1
fi

# Dizinleri oluştur
mkdir -p "$OUT_DIR"

# Resume durumu kontrol et
check_resume_status

# Zaman takibi
START_TIME=$(date +%s)
START_TIME_FMT=$(date '+%Y-%m-%d %H:%M:%S')
TOTAL_FILES=$(wc -l < "$ACC_LIST")

echo "============================================================"
echo "🚀 GNU PARALLEL İLE FASTQ DÖNÜŞTÜRME İŞLEMİ"
echo "============================================================"
echo "⏰ Başlangıç: $START_TIME_FMT"
echo "📁 Toplam Dosya: $TOTAL_FILES"
echo "🔄 Paralel İşlem: $MAX_PARALLEL CPU çekirdeği"
echo "🔁 Maksimum Retry: $MAX_RETRIES"
echo "📂 Çıktı Dizini: $OUT_DIR"
echo "📝 Detaylı Log: $PARALLEL_LOG"
echo "💾 Minimum Dosya Boyutu: ${MIN_FILE_SIZE} bytes"
echo "============================================================"

# --- GNU PARALLEL İLE ANA İŞLEME ---
cat "$ACC_LIST" | parallel \
    --jobs "$MAX_PARALLEL" \
    --retries "$MAX_RETRIES" \
    --joblog "$PARALLEL_LOG" \
    --resume \
    --bar \
    --halt-on-error 0 \
    --timeout 2000 \
    'process_srr {} {#} '"$TOTAL_FILES"

# Parallel'in exit kodunu kaydet
PARALLEL_EXIT=$?

# --- İŞLEM SONU RAPORLAMA ---
END_TIME=$(date +%s)
TOTAL_SECONDS=$((END_TIME - START_TIME))
TOTAL_TIME=$(date -u -d @${TOTAL_SECONDS} +'%H:%M:%S' 2>/dev/null || echo "${TOTAL_SECONDS}s")

# İstatistikleri joblog'dan hesapla
if [ -f "$PARALLEL_LOG" ]; then
    SUCCESS_COUNT=$(awk 'NR > 1 && $7 == 0' "$PARALLEL_LOG" | wc -l)
    FAIL_COUNT=$(awk 'NR > 1 && $7 != 0' "$PARALLEL_LOG" | wc -l)
    TOTAL_PROCESSED=$((SUCCESS_COUNT + FAIL_COUNT))
    
    # Hatalı SRR'ları log dosyasına yaz
    awk 'NR > 1 && $7 != 0 {print $9}' "$PARALLEL_LOG" > "$LOG_FILE"
else
    SUCCESS_COUNT=0
    FAIL_COUNT=0
    TOTAL_PROCESSED=0
fi

# Çıktı dizinindeki gerçek dosya sayısını say
ACTUAL_FILES=$(find "$OUT_DIR" -name "*.fastq.gz" | wc -l)
TOTAL_SIZE=$(du -sh "$OUT_DIR" 2>/dev/null | cut -f1 || echo "N/A")

echo -e "\n\n============================================================"
echo "📊 İŞLEM RAPORU"
echo "============================================================"
echo "⏰ Toplam Süre: $TOTAL_TIME"
echo "📁 İşlenen Dosya: $TOTAL_PROCESSED / $TOTAL_FILES"
echo "✅ Başarılı: $SUCCESS_COUNT"
echo "❌ Başarısız: $FAIL_COUNT"
echo "📈 Başarı Oranı: $(( SUCCESS_COUNT * 100 / TOTAL_FILES ))%"
echo "💾 Çıktı Dosyaları: $ACTUAL_FILES adet ($TOTAL_SIZE)"

if [ "$TOTAL_SECONDS" -gt 0 ]; then
    echo "⚡ Ortalama Hız: $(( TOTAL_PROCESSED * 3600 / TOTAL_SECONDS )) dosya/saat"
fi

if [ "$FAIL_COUNT" -gt 0 ] && [ -f "$LOG_FILE" ]; then
    echo "------------------------------------------------------------"
    echo "💥 Başarısız SRR ID'leri ($FAIL_COUNT adet):"
    head -5 "$LOG_FILE" | sed 's/^/   /'
    [ "$FAIL_COUNT" -gt 5 ] && echo "   ... ve $(( FAIL_COUNT - 5 )) adet daha"
    echo ""
    echo "📝 Tam liste: $LOG_FILE"
    echo "🔄 Yalnızca başarısızları yeniden çalıştırmak için:"
    echo "   cat $LOG_FILE | parallel [aynı parametreler]"
fi

echo "============================================================"

# Çıkış kodu belirleme
if [ "$FAIL_COUNT" -eq 0 ] && [ "$PARALLEL_EXIT" -eq 0 ]; then
    echo "🎉 Tüm işlemler başarıyla tamamlandı!"
    exit 0
elif [ "$SUCCESS_COUNT" -gt 0 ]; then
    echo "⚠️ İşlem kısmen başarılı. Başarısızlar için log dosyasını kontrol edin."
    exit 1
else
    echo "❌ İşlem başarısız. Lütfen konfigürasyonu ve log dosyalarını kontrol edin."
    exit 2
fi