"""
============================================================================
Genome Analysis Pipeline - Main Snakefile
============================================================================
A modular Snakemake pipeline for bacterial and viral genome analysis.

Modules:
    1. Quality Control (QC)
    2. Assembly
    3. Identification
    4. Annotation
    5. Specialized Analysis (AMR, Virulence, MGE, Prophages)
    6. Anomaly Detection
    7. (Reserved for future use)
    8. Report Generation

Usage:
    snakemake --configfile config/config.yaml --cores 8 --use-conda

    # With custom sample:
    snakemake --configfile config/config.yaml \
        --config sample_id=my_sample \
                 reads.short.r1=reads_R1.fastq.gz \
                 reads.short.r2=reads_R2.fastq.gz \
        --cores 8 --use-conda
"""

import os
import sys
from pathlib import Path

# ============================================================================
# Configuration
# ============================================================================

# Load main config
configfile: "config/config.yaml"

# Load database paths if available
if os.path.exists("config/databases.yaml"):
    configfile: "config/databases.yaml"

# Include common functions
include: "workflow/rules/common.smk"

# ============================================================================
# Global Variables
# ============================================================================

SAMPLE_ID = get_sample_id(config)
OUTPUT_DIR = get_sample_output(config)
THREADS = get_threads(config)

# Detect read mode (only required if running upstream modules)
def safe_detect_read_mode():
    """Detect read mode, return 'unknown' if reads not configured"""
    try:
        return detect_read_mode(config)
    except ValueError:
        # If only running report module, reads aren't needed
        if is_module_enabled(config, "qc") or is_module_enabled(config, "assembly"):
            raise
        return "unknown"

READ_MODE = safe_detect_read_mode()

# ============================================================================
# Include Module Rules
# ============================================================================

include: "workflow/rules/qc.smk"
include: "workflow/rules/assembly.smk"
include: "workflow/rules/identification.smk"
include: "workflow/rules/annotation.smk"
include: "workflow/rules/amr.smk"
include: "workflow/rules/virulence.smk"
include: "workflow/rules/mobile_elements.smk"
include: "workflow/rules/prophages.smk"
include: "workflow/rules/anomaly.smk"
include: "workflow/rules/report.smk"

# ============================================================================
# Target Rules
# ============================================================================

def get_final_outputs():
    """Determine final outputs based on enabled modules"""
    outputs = []
    
    # Always include logs
    outputs.append(f"{OUTPUT_DIR}/00_logs/pipeline_complete.log")
    
    # QC outputs
    if is_module_enabled(config, "qc"):
        outputs.append(f"{OUTPUT_DIR}/01_qc/multiqc_report.html")
    
    # Assembly outputs
    if is_module_enabled(config, "assembly"):
        outputs.append(f"{OUTPUT_DIR}/02_assembly/contigs.fasta")
        outputs.append(f"{OUTPUT_DIR}/02_assembly/quast/report.html")
    
    # Identification outputs
    if is_module_enabled(config, "identification"):
        outputs.append(f"{OUTPUT_DIR}/03_identification/taxonomy.json")
        outputs.append(f"{OUTPUT_DIR}/03_identification/organism_type.txt")
    
    # Annotation outputs
    if is_module_enabled(config, "annotation"):
        outputs.append(f"{OUTPUT_DIR}/04_annotation/annotation_complete.flag")
    
    # Specialized outputs
    if is_module_enabled(config, "amr"):
        outputs.append(f"{OUTPUT_DIR}/05_specialized/amr/amrfinder_results.tsv")
    
    if is_module_enabled(config, "virulence"):
        outputs.append(f"{OUTPUT_DIR}/05_specialized/virulence/vfdb_results.tsv")
    
    if is_module_enabled(config, "mobile_elements"):
        outputs.append(f"{OUTPUT_DIR}/05_specialized/mge/mob_suite_complete.flag")
    
    if is_module_enabled(config, "prophages"):
        outputs.append(f"{OUTPUT_DIR}/05_specialized/prophages/phispy_complete.flag")
    
    # Anomaly detection
    if is_module_enabled(config, "anomaly"):
        outputs.append(f"{OUTPUT_DIR}/06_anomaly/anomaly_summary.json")
    
    # Report
    if is_module_enabled(config, "report"):
        outputs.append(f"{OUTPUT_DIR}/08_report/genome_report.html")
    
    return outputs


rule all:
    """Main target rule - runs entire pipeline"""
    input:
        get_final_outputs()


rule clean:
    """Remove all output files for this sample"""
    shell:
        f"rm -rf {OUTPUT_DIR}"


# ============================================================================
# Utility Rules
# ============================================================================

rule create_output_dirs:
    """Create output directory structure"""
    output:
        directory(f"{OUTPUT_DIR}/00_logs"),
        directory(f"{OUTPUT_DIR}/01_qc"),
        directory(f"{OUTPUT_DIR}/02_assembly"),
        directory(f"{OUTPUT_DIR}/03_identification"),
        directory(f"{OUTPUT_DIR}/04_annotation"),
        directory(f"{OUTPUT_DIR}/05_specialized"),
        directory(f"{OUTPUT_DIR}/06_anomaly"),
        directory(f"{OUTPUT_DIR}/08_report"),
    shell:
        """
        mkdir -p {output}
        mkdir -p {OUTPUT_DIR}/05_specialized/{{amr,virulence,mge,prophages}}
        mkdir -p {OUTPUT_DIR}/06_anomaly/{{reference,alignment,gene_taxonomy}}
        """


rule log_pipeline_start:
    """Log pipeline start information"""
    output:
        f"{OUTPUT_DIR}/00_logs/pipeline_start.log"
    run:
        import datetime
        with open(output[0], "w") as f:
            f.write(f"Pipeline started: {datetime.datetime.now()}\n")
            f.write(f"Sample ID: {SAMPLE_ID}\n")
            f.write(f"Read mode: {READ_MODE}\n")
            f.write(f"Threads: {THREADS}\n")
            f.write(f"Config: {config}\n")


rule log_pipeline_complete:
    """Log pipeline completion"""
    input:
        # Depend on report as final output
        report = f"{OUTPUT_DIR}/08_report/genome_report.html" if is_module_enabled(config, "report") else [],
    output:
        f"{OUTPUT_DIR}/00_logs/pipeline_complete.log"
    run:
        import datetime
        with open(output[0], "w") as f:
            f.write(f"Pipeline completed: {datetime.datetime.now()}\n")
            f.write(f"Sample ID: {SAMPLE_ID}\n")
            f.write("All modules finished successfully.\n")


# ============================================================================
# Pipeline Information
# ============================================================================

onstart:
    print("=" * 70)
    print("GENOME ANALYSIS PIPELINE")
    print("=" * 70)
    print(f"Sample ID: {SAMPLE_ID}")
    print(f"Read mode: {READ_MODE}")
    print(f"Output: {OUTPUT_DIR}")
    print(f"Threads: {THREADS}")
    print("=" * 70)


onsuccess:
    print("=" * 70)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print("=" * 70)
    print(f"Results in: {OUTPUT_DIR}")
    if is_module_enabled(config, "report"):
        print(f"Report: {OUTPUT_DIR}/08_report/genome_report.html")
    print("=" * 70)


onerror:
    print("=" * 70)
    print("PIPELINE FAILED")
    print("=" * 70)
    print(f"Check logs in: {OUTPUT_DIR}/00_logs/")
    print("=" * 70)

