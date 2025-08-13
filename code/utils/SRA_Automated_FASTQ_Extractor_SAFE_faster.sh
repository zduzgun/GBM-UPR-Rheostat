#!/bin/bash

# ====================================================================
# GNU PARALLEL + FASTERQ-DUMP + PIGZ İLE GÜÇLENDİRİLMİŞ
# SRA->FASTQ DÖNÜŞTÜRME - ULTRA HIZLI VE AKILLI VERSİYON
#
# YENİ ÖZELLİK: Otomatik Liste Güncelleme
# Betik her çalıştığında veya durdurulduğunda, başarılı dosyaları
# listeden otomatik olarak çıkarır ve bir sonraki çalışmada
# sadece kalanlardan devam eder.
# ====================================================================

# --- YAPILANDIRMA DEĞİŞKENLERİ ---
OUT_DIR="fastq_files"
ACC_LIST="resmi_liste_875.txt" # Ana çalışma listesi, bu dosya güncellenecektir.
LOG_FILE="hatali_dosyalar.log"
PARALLEL_LOG="parallel_job.log"
TEMP_DIR="/tmp/sra_temp_$$"  # Process-specific temp dizini

# Orijinal listenin yedeği için dosya adı
MASTER_LIST_BACKUP="${ACC_LIST}.orig"

# Paralel işlem sayısı
MAX_PARALLEL=$(nproc)
MAX_RETRIES=3
MIN_FILE_SIZE=1024  # 1KB minimum dosya boyutu

# PIGZ için thread sayısı (CPU çekirdeği başına)
PIGZ_THREADS=2

# --- YARDIMCI FONKSİYONLAR ---

# YENİ FONKSİYON: Başarılı dosyaları analiz ederek ACC_LIST'i günceller
update_acc_list() {
    echo -e "\n🔄 Aksesyon listesi kontrol ediliyor ve güncelleniyor..."

    # Orijinal liste yedeği mevcut değilse işlem yapma (güvenlik)
    if [ ! -f "$MASTER_LIST_BACKUP" ]; then
        echo "⚠️  Orijinal liste yedeği ($MASTER_LIST_BACKUP) bulunamadı. Liste güncellenmedi."
        return
    fi
    
    # Çıktı dizini yoksa veya boşsa, işlem yapmaya gerek yok.
    if [ ! -d "$OUT_DIR" ] || [ -z "$(ls -A "$OUT_DIR" 2>/dev/null)" ]; then
        echo "✅ Çıktı dizininde işlenmiş dosya bulunamadı. Liste aynı kaldı."
        return
    fi
    
    local success_file="${TEMP_DIR}/success_srr.txt"
    local temp_acc_list="${TEMP_DIR}/new_acc_list.txt"

    # 1. Başarılı şekilde oluşturulmuş dosyaların SRR ID'lerini al
    ls -1 "$OUT_DIR" 2>/dev/null | sed -E 's/(_[12])?\.fastq\.gz$//' | sort | uniq > "$success_file"

    # 2. Mevcut aksesyon listesinden başarılı olanları çıkar
    if [ -s "$success_file" ]; then
        grep -vFf "$success_file" "$ACC_LIST" > "$temp_acc_list"
        
        local original_count=$(wc -l < "$ACC_LIST")
        local new_count=$(wc -l < "$temp_acc_list")
        
        # 3. Ana listeyi güncelle
        mv "$temp_acc_list" "$ACC_LIST"
        
        echo "✅ Liste güncellendi. $original_count dosyadan $new_count dosya kaldı."
    else
        echo "✅ Yeni başarılı dosya bulunamadı. Liste değiştirilmedi."
    fi
}

# Temizlik fonksiyonu (geliştirilmiş ve liste güncelleme eklenmiş)
cleanup() {
    echo -e "\n\n🛑 İşlem durduruldu veya tamamlandı. Temizlik ve güncelleme yapılıyor..."
    
    # Aktif parallel işlemlerini durdur
    pkill -f "parallel.*process_srr" 2>/dev/null
    pkill -f "fasterq-dump" 2>/dev/null
    pkill -f "pigz" 2>/dev/null
    
    # --- YENİ EKLENEN KISIM ---
    # Aksesyon listesini güncelle, böylece bir sonraki çalıştırmada sadece eksikler kalır
    update_acc_list
    # --------------------------
    
    # Temporary dosyaları temizle
    echo "🧼 Geçici dosyalar temizleniyor..."
    rm -rf "$TEMP_DIR" 2>/dev/null
    find "$OUT_DIR" -name "*.tmp" -delete 2>/dev/null
    find "$OUT_DIR" -name "*.fastq" -delete 2>/dev/null  # Sıkıştırılmamış fastq'lar
    
    echo "✅ Temizlik tamamlandı. Betiği yeniden çalıştırarak kalanlardan devam edebilirsiniz."
}

# Sinyal yakalama (EXIT: normal bitişte de çalışır, SIGINT: Ctrl+C'de çalışır)
trap cleanup EXIT SIGINT SIGTERM

# Resume kontrolü
check_resume_status() {
    if [ -f "$PARALLEL_LOG" ]; then
        local completed=$(awk 'NR > 1 && $7 == 0' "$PARALLEL_LOG" | wc -l)
        local failed=$(awk 'NR > 1 && $7 != 0' "$PARALLEL_LOG" | wc -l)
        local total=$(wc -l < "$ACC_LIST")
        local remaining=$((total - completed - failed))
        
        if [ "$remaining" -gt 0 ]; then
            echo "🔄 Önceki parallel.log tespit edildi:"
            echo "   ✅ Tamamlanan: $completed"
            echo "   ❌ Başarısız: $failed"
            echo "   (Bu bilgiler sadece referans içindir, betik kalanları kendi listesinden okuyacaktır)"
        fi
    fi
    return 1
}

# Çıktı dosyası doğrulama
verify_output() {
    local srr=$1
    local verified=0
    
    for suffix in ".fastq.gz" "_1.fastq.gz" "_2.fastq.gz"; do
        local file="${OUT_DIR}/${srr}${suffix}"
        if [ -f "$file" ]; then
            local size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null || echo 0)
            if [ "$size" -gt "$MIN_FILE_SIZE" ]; then
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

# Ana işleme fonksiyonu (fasterq-dump + pigz optimize)
process_srr() {
    local srr=$1
    local job_num=$2
    local total_jobs=$3
    
    local job_temp_dir="${TEMP_DIR}/${srr}_$$"
    mkdir -p "$job_temp_dir"
    
    if [ -f "${OUT_DIR}/${srr}_1.fastq.gz" ] || [ -f "${OUT_DIR}/${srr}.fastq.gz" ]; then
        echo "($job_num/$total_jobs) ⏭️ $srr zaten mevcut, atlanıyor."
        rm -rf "$job_temp_dir"
        return 0
    fi
    
    echo "($job_num/$total_jobs) 🔄 $srr işleniyor... (fasterq-dump + pigz)"
    
    if timeout 1200 fasterq-dump "$srr" \
        --outdir "$job_temp_dir" \
        --temp "$job_temp_dir" \
        --split-3 \
        --skip-technical \
        --progress 2>/dev/null; then
        
        echo "($job_num/$total_jobs) 📦 $srr FASTQ çıkarıldı, sıkıştırılıyor..."
        
        local compress_success=true
        
        for fastq_file in "$job_temp_dir"/*.fastq; do
            if [ -f "$fastq_file" ]; then
                local base_name=$(basename "$fastq_file" .fastq)
                local output_file="${OUT_DIR}/${base_name}.fastq.gz"
                
                if ! pigz -p "$PIGZ_THREADS" -c "$fastq_file" > "$output_file"; then
                    echo "⚠️ $fastq_file sıkıştırma hatası"
                    compress_success=false
                    rm -f "$output_file"
                fi
                rm -f "$fastq_file"
            fi
        done
        
        rm -rf "$job_temp_dir"
        
        if [ "$compress_success" = true ] && verify_output "$srr"; then
            echo "($job_num/$total_jobs) ✅ $srr başarıyla dönüştürüldü ve sıkıştırıldı."
            return 0
        else
            echo "($job_num/$total_jobs) ❌ $srr sıkıştırma veya doğrulama başarısız."
            rm -f "${OUT_DIR}/${srr}"*.fastq.gz 2>/dev/null
            return 1
        fi
    else
        echo "($job_num/$total_jobs) ❌ $srr fasterq-dump hatası veya timeout."
        rm -rf "$job_temp_dir"
        rm -f "${OUT_DIR}/${srr}"*.fastq.gz "${OUT_DIR}/${srr}"*.tmp 2>/dev/null
        return 1
    fi
}

# Fonksiyonları GNU Parallel için export et
export -f process_srr verify_output update_acc_list
export OUT_DIR MIN_FILE_SIZE TEMP_DIR PIGZ_THREADS ACC_LIST MASTER_LIST_BACKUP

# ====================================================================
# ANA PROGRAM
# ====================================================================

# Başlangıç kontrolleri
if [ ! -f "$ACC_LIST" ]; then
    if [ -f "$MASTER_LIST_BACKUP" ]; then
        echo "⚠️  $ACC_LIST bulunamadı, yedekten ($MASTER_LIST_BACKUP) geri yükleniyor..."
        cp "$MASTER_LIST_BACKUP" "$ACC_LIST"
    else
        echo "❌ HATA: $ACC_LIST dosyası ve yedeği bulunamadı!"
        exit 1
    fi
fi

# YENİ: Orijinal listenin yedeğini oluştur (eğer yoksa)
if [ ! -f "$MASTER_LIST_BACKUP" ]; then
    echo "✨ İlk çalıştırma tespit edildi. Orijinal listenin yedeği oluşturuluyor: $MASTER_LIST_BACKUP"
    cp "$ACC_LIST" "$MASTER_LIST_BACKUP"
fi

# Gerekli araçları kontrol et
missing_tools=()
for tool in fasterq-dump pigz parallel; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        missing_tools+=("$tool")
    fi
done

if [ ${#missing_tools[@]} -gt 0 ]; then
    echo "❌ Hata: Aşağıdaki araçlar bulunamadı:"
    printf '   - %s\n' "${missing_tools[@]}"
    echo "Lütfen SRA Toolkit, pigz ve GNU Parallel kurun."
    exit 1
fi

# Dizinleri oluştur
mkdir -p "$OUT_DIR"
mkdir -p "$TEMP_DIR"

# Liste boş mu kontrol et
if ! [ -s "$ACC_LIST" ]; then
    echo -e "\n🎉 Tebrikler! Aksesyon listesi ($ACC_LIST) boş."
    echo "Tüm dosyalar başarıyla işlenmiş görünüyor."
    echo "Yeniden başlamak için orijinal yedeği geri yükleyebilirsiniz:"
    echo "cp '$MASTER_LIST_BACKUP' '$ACC_LIST'"
    exit 0
fi

# `parallel` log'unu temiz bir başlangıç için sıfırla (isteğe bağlı ama önerilir)
# Bu satırı aktif ederseniz her seferinde parallel'in log'u silinir.
# Otomatik liste güncelleme olduğu için bu artık daha az önemli.
# > "$PARALLEL_LOG"

# check_resume_status # Bu fonksiyon artık sadece bilgilendirme amaçlı

# Zaman takibi
START_TIME=$(date +%s)
START_TIME_FMT=$(date '+%Y-%m-%d %H:%M:%S')
TOTAL_FILES=$(wc -l < "$ACC_LIST")

echo "============================================================"
echo "🚀 AKILLI SRA DÖNÜŞTÜRÜCÜ (Kaldığı Yerden Devam Eder)"
echo "============================================================"
echo "⏰ Başlangıç: $START_TIME_FMT"
echo "📁 İşlenecek Dosya: $TOTAL_FILES"
echo "🔄 Paralel İşlem: $MAX_PARALLEL CPU çekirdeği"
echo "🗜️ PIGZ Thread/İş: $PIGZ_THREADS"
echo "📂 Çıktı Dizini: $OUT_DIR"
echo "📝 Detaylı Log: $PARALLEL_LOG"
echo "============================================================"

# SRA toolkit konfigürasyonu
vdb-config --set /repository/user/cache-enabled=true 2>/dev/null
vdb-config --set /repository/user/cache-size=10G 2>/dev/null

# --- GNU PARALLEL İLE ANA İŞLEME ---
cat "$ACC_LIST" | parallel \
    --jobs "$MAX_PARALLEL" \
    --retries "$MAX_RETRIES" \
    --joblog "$PARALLEL_LOG" \
    --bar \
    --halt-on-error 0 \
    --timeout 1800 \
    --memfree 1G \
    'process_srr {} {#} '"$TOTAL_FILES"

# Parallel'in exit kodunu kaydet
PARALLEL_EXIT=$?

# --- İŞLEM SONU RAPORLAMA ---
# Not: `trap` komutu nedeniyle bu bölümden sonra `cleanup` fonksiyonu çalışacaktır.
END_TIME=$(date +%s)
TOTAL_SECONDS=$((END_TIME - START_TIME))
TOTAL_TIME=$(date -u -d @${TOTAL_SECONDS} +'%H:%M:%S' 2>/dev/null || echo "${TOTAL_SECONDS}s")

if [ -f "$PARALLEL_LOG" ]; then
    SUCCESS_COUNT=$(awk 'NR > 1 && $7 == 0' "$PARALLEL_LOG" | wc -l)
else
    SUCCESS_COUNT=0
fi

echo -e "\n\n============================================================"
echo "📊 İŞLEM TAMAMLANDI RAPORU"
echo "============================================================"
echo "⏰ Toplam Süre: $TOTAL_TIME"
echo "✅ Bu oturumda başarıyla işlenen dosya: $SUCCESS_COUNT"

# Çıkış kodu belirleme
if [ "$PARALLEL_EXIT" -eq 0 ]; then
    echo "🎉 Bu oturumdaki tüm işlemler başarıyla tamamlandı."
    exit 0
else
    echo "⚠️  Bu oturumda bazı hatalar oluştu. Liste güncellendi, kalanlar için yeniden çalıştırın."
    exit 1
fi