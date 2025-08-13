#!/usr/bin/env python3
"""
Enterprise Kallisto Pipeline - Configuration Examples
====================================================

Bu dosya farklı kullanım senaryoları için örnek konfigürasyonlar içerir.
Kendi kullanımınız için uygun olanı seçip düzenleyebilirsiniz.
"""

from pathlib import Path
from run_kallisto import PipelineConfig, TechnologyType, EnterpriseKallistoPipeline

# =============================================================================
# TEMEL KONFİGÜRASYONLAR
# =============================================================================

# Smart-seq2 - Küçük veri seti için
SMARTSEQ2_SMALL = PipelineConfig(
    fastq_dir=Path("fastq_files"),
    output_dir=Path("results/smartseq2_small"),
    index_path=Path("reference/human_transcriptome.idx"),
    technology=TechnologyType.SMART_SEQ2,
    max_workers=4,
    threads_per_sample=4,
    timeout_minutes=30,
    memory_limit_gb=8,
    max_retries=3,
    health_check_interval=30,
    log_level="INFO"
)

# Smart-seq2 - Orta boy veri seti için
SMARTSEQ2_MEDIUM = PipelineConfig(
    fastq_dir=Path("fastq_files"),
    output_dir=Path("results/smartseq2_medium"),
    index_path=Path("reference/human_transcriptome.idx"),
    technology=TechnologyType.SMART_SEQ2,
    max_workers=16,
    threads_per_sample=4,
    timeout_minutes=60,
    memory_limit_gb=32,
    max_retries=5,
    retry_delay_seconds=60,
    health_check_interval=30,
    checkpoint_interval=300,
    log_level="INFO"
)

# Smart-seq2 - Production ortamı için
SMARTSEQ2_PRODUCTION = PipelineConfig(
    fastq_dir=Path("/data/production/fastq"),
    output_dir=Path("/results/production/smartseq2"),
    index_path=Path("/reference/production/human_transcriptome.idx"),
    technology=TechnologyType.SMART_SEQ2,
    max_workers=64,
    threads_per_sample=2,
    timeout_minutes=120,
    memory_limit_gb=128,
    max_retries=10,
    retry_delay_seconds=120,
    health_check_interval=60,
    checkpoint_interval=300,
    enable_monitoring=True,
    auto_scaling=True,
    log_level="INFO",
    log_rotation_mb=500,
    max_log_files=20
)

# 10x Chromium v3 - Standart konfigürasyon
TENX_V3_STANDARD = PipelineConfig(
    fastq_dir=Path("10x_fastq"),
    output_dir=Path("results/10x_v3"),
    index_path=Path("reference/human_transcriptome.idx"),
    technology=TechnologyType.TENX_V3,
    max_workers=16,
    threads_per_sample=4,
    timeout_minutes=90,
    memory_limit_gb=32,
    run_bustools=True,
    generate_count_matrix=True,
    max_retries=5,
    health_check_interval=30,
    log_level="INFO"
)

# Drop-seq - Araştırma ortamı için
DROPSEQ_RESEARCH = PipelineConfig(
    fastq_dir=Path("dropseq_data"),
    output_dir=Path("results/dropseq"),
    index_path=Path("reference/mouse_transcriptome.idx"),
    technology=TechnologyType.DROP_SEQ,
    max_workers=8,
    threads_per_sample=6,
    timeout_minutes=120,
    memory_limit_gb=24,
    run_bustools=True,
    max_retries=3,
    log_level="DEBUG"
)

# =============================================================================
# ÖZEL KULLANIM SENARYOLARI
# =============================================================================

# Yüksek bellek optimizasyonu - Bellek kısıtlı sistemler için
HIGH_MEMORY_OPTIMIZATION = PipelineConfig(
    fastq_dir=Path("fastq_files"),
    output_dir=Path("results/memory_optimized"),
    index_path=Path("reference/human_transcriptome.idx"),
    technology=TechnologyType.SMART_SEQ2,
    max_workers=4,  # Düşük işçi sayısı
    threads_per_sample=8,  # Yüksek thread sayısı
    timeout_minutes=180,  # Uzun timeout
    memory_limit_gb=16,
    max_retries=5,
    auto_scaling=False,  # Manuel kontrol
    cleanup_temp=True,
    log_level="WARNING"  # Düşük log seviyesi
)

# Hızlı işleme - Yüksek performans sistemleri için
FAST_PROCESSING = PipelineConfig(
    fastq_dir=Path("fastq_files"),
    output_dir=Path("results/fast_processing"),
    index_path=Path("reference/human_transcriptome.idx"),
    technology=TechnologyType.SMART_SEQ2,
    max_workers=32,  # Yüksek işçi sayısı
    threads_per_sample=2,  # Düşük thread sayısı
    timeout_minutes=30,  # Kısa timeout
    memory_limit_gb=64,
    max_retries=3,
    retry_delay_seconds=30,
    health_check_interval=15,  # Sık health check
    checkpoint_interval=180,  # Sık checkpoint
    auto_scaling=True,
    log_level="INFO"
)

# Debug ve geliştirme ortamı
DEBUG_DEVELOPMENT = PipelineConfig(
    fastq_dir=Path("test_data/fastq"),
    output_dir=Path("test_results/debug"),
    index_path=Path("test_data/test_index.idx"),
    technology=TechnologyType.SMART_SEQ2,
    max_workers=2,
    threads_per_sample=2,
    timeout_minutes=15,
    memory_limit_gb=4,
    max_retries=1,
    dry_run=False,
    validate_inputs=True,
    cleanup_temp=False,  # Temp dosyaları tut
    log_level="DEBUG",
    health_check_interval=10
)

# =============================================================================
# KULLANIM ÖRNEKLERİ
# =============================================================================

def run_smartseq2_small():
    """Küçük Smart-seq2 veri seti için pipeline çalıştır"""
    pipeline = EnterpriseKallistoPipeline(SMARTSEQ2_SMALL)
    report = pipeline.run_enterprise_pipeline()
    pipeline.print_enterprise_summary(report)
    return report

def run_production_pipeline():
    """Production ortamında pipeline çalıştır"""
    pipeline = EnterpriseKallistoPipeline(SMARTSEQ2_PRODUCTION)
    report = pipeline.run_enterprise_pipeline()
    pipeline.print_enterprise_summary(report)
    return report

def run_10x_analysis():
    """10x Chromium v3 analizi çalıştır"""
    pipeline = EnterpriseKallistoPipeline(TENX_V3_STANDARD)
    report = pipeline.run_enterprise_pipeline()
    pipeline.print_enterprise_summary(report)
    return report

def dry_run_test():
    """Dry run testi yap"""
    config = SMARTSEQ2_SMALL
    config.dry_run = True
    
    pipeline = EnterpriseKallistoPipeline(config)
    report = pipeline.run_enterprise_pipeline()
    print("Dry run tamamlandı!")
    return report

# =============================================================================
# KONFİGÜRASYON YARDIMCI FONKSİYONLARI
# =============================================================================

def get_config_for_system_specs(cpu_cores: int, memory_gb: int, sample_count: int) -> PipelineConfig:
    """
    Sistem özelliklerine göre optimal konfigürasyon oluştur
    
    Args:
        cpu_cores: CPU çekirdek sayısı
        memory_gb: Sistem belleği (GB)
        sample_count: Örnek sayısı
    
    Returns:
        Optimize edilmiş PipelineConfig
    """
    # Temel hesaplamalar
    max_workers = min(cpu_cores // 2, memory_gb // 4)  # Bellek ve CPU'ya göre
    threads_per_sample = max(2, cpu_cores // max_workers)
    memory_limit = memory_gb * 0.8  # Sistem belleğinin %80'i
    
    # Örnek sayısına göre ayarlamalar
    if sample_count > 1000:
        # Büyük veri seti - güvenilirlik odaklı
        max_retries = 10
        timeout_minutes = 180
        checkpoint_interval = 300
    elif sample_count > 100:
        # Orta boy veri seti - standart ayarlar
        max_retries = 5
        timeout_minutes = 90
        checkpoint_interval = 600
    else:
        # Küçük veri seti - hızlı işleme
        max_retries = 3
        timeout_minutes = 60
        checkpoint_interval = 900
    
    return PipelineConfig(
        fastq_dir=Path("fastq_files"),
        output_dir=Path("results/auto_optimized"),
        index_path=Path("reference/human_transcriptome.idx"),
        technology=TechnologyType.SMART_SEQ2,
        max_workers=max_workers,
        threads_per_sample=threads_per_sample,
        timeout_minutes=timeout_minutes,
        memory_limit_gb=memory_limit,
        max_retries=max_retries,
        checkpoint_interval=checkpoint_interval,
        auto_scaling=True,
        log_level="INFO"
    )

def create_custom_config(**kwargs) -> PipelineConfig:
    """
    Özel parametrelerle konfigürasyon oluştur
    
    Args:
        **kwargs: PipelineConfig parametreleri
        
    Returns:
        Özelleştirilmiş PipelineConfig
    """
    # Varsayılan değerler
    defaults = {
        'fastq_dir': Path("fastq_files"),
        'output_dir': Path("results/custom"),
        'index_path': Path("reference/human_transcriptome.idx"),
        'technology': TechnologyType.SMART_SEQ2,
        'max_workers': 8,
        'threads_per_sample': 4,
        'timeout_minutes': 60,
        'memory_limit_gb': 16,
        'max_retries': 3,
        'log_level': "INFO"
    }
    
    # Kullanıcı parametreleri ile güncelle
    defaults.update(kwargs)
    
    return PipelineConfig(**defaults)

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Enterprise Kallisto Pipeline - Konfigürasyon Örnekleri"
    )
    parser.add_argument(
        "--config",
        choices=["small", "medium", "production", "10x", "dropseq", "debug"],
        default="small",
        help="Kullanılacak konfigürasyon"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Dry run modu"
    )
    
    args = parser.parse_args()
    
    # Konfigürasyon seçimi
    config_map = {
        "small": SMARTSEQ2_SMALL,
        "medium": SMARTSEQ2_MEDIUM,
        "production": SMARTSEQ2_PRODUCTION,
        "10x": TENX_V3_STANDARD,
        "dropseq": DROPSEQ_RESEARCH,
        "debug": DEBUG_DEVELOPMENT
    }
    
    config = config_map[args.config]
    
    if args.dry_run:
        config.dry_run = True
    
    print(f"🚀 {args.config.upper()} konfigürasyonu ile pipeline başlatılıyor...")
    print(f"📋 Konfigürasyon detayları:")
    print(f"   • FASTQ Directory: {config.fastq_dir}")
    print(f"   • Output Directory: {config.output_dir}")
    print(f"   • Technology: {config.technology.value}")
    print(f"   • Workers: {config.max_workers}")
    print(f"   • Memory Limit: {config.memory_limit_gb} GB")
    print(f"   • Max Retries: {config.max_retries}")
    print(f"   • Dry Run: {config.dry_run}")
    
    # Pipeline çalıştır
    pipeline = EnterpriseKallistoPipeline(config)
    report = pipeline.run_enterprise_pipeline()
    pipeline.print_enterprise_summary(report) 