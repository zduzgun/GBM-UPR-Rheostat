#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
from run_kallisto import PipelineConfig, TechnologyType, EnterpriseKallistoPipeline

def main():
    config = PipelineConfig(
        fastq_dir=Path("fastq_files"),
        output_dir=Path("test_output"),
        index_path=Path("referans/insan_transkriptom.idx"),
        technology=TechnologyType.SMART_SEQ2,
        max_workers=2,
        threads_per_sample=2,
        timeout_minutes=30,
        memory_limit_gb=16,
        max_retries=1,
        dry_run=True,
        enable_monitoring=False,
        log_level="INFO"
    )
    
    print("Starting Enterprise Kallisto Pipeline Test")
    print(f"Configuration: {config.technology.value}")
    print(f"Workers: {config.max_workers}")
    print(f"Dry run: {config.dry_run}")
    
    try:
        pipeline = EnterpriseKallistoPipeline(config)
        report = pipeline.run_enterprise_pipeline()
        pipeline.print_enterprise_summary(report)
        print("Test completed successfully!")
        
    except Exception as e:
        print(f"Test failed: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 