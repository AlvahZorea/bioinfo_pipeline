# Comprehensive Bacterial and Viral Genome Analysis Pipeline

## Executive Summary

This document outlines the design and preparation for a modular bioinformatics pipeline for deep characterization of bacterial and viral genomes. The pipeline will:

- Accept Illumina short reads, Oxford Nanopore long reads, or hybrid input
- Handle pure cultures and mixed samples
- Support bacteria and all viruses (DNA, RNA, bacteriophages)
- Run in fully offline (air-gapped) environments
- Generate comprehensive HTML reports with interactive visualizations
- Be highly modular and extensible using Snakemake

## Pipeline Architecture

```
INPUT                    CORE PROCESSING                         OUTPUT
─────                    ───────────────                         ──────

┌─────────┐     ┌────────────────────────────────────┐     ┌─────────────┐
│ FASTQ   │────▶│  MODULE 1: Quality Control         │────▶│ QC Report   │
│ (Short) │     │  - FastQC/MultiQC                  │     │             │
│         │     │  - Trimming/Filtering              │     │             │
└─────────┘     │  - Contamination screening         │     └─────────────┘
                └────────────────────────────────────┘
┌─────────┐                      │
│ FASTQ   │                      ▼
│ (Long)  │     ┌────────────────────────────────────┐     ┌─────────────┐
│         │────▶│  MODULE 2: Assembly                │────▶│ Contigs     │
└─────────┘     │  - Short-read assembly             │     │ FASTA       │
                │  - Long-read assembly              │     │             │
┌─────────┐     │  - Hybrid assembly                 │     └─────────────┘
│ BOTH    │────▶│  - Metagenome binning (if mixed)   │
└─────────┘     └────────────────────────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐     ┌─────────────┐
                │  MODULE 3: Identification          │────▶│ Taxonomy    │
                │  - Taxonomic classification        │     │ Report      │
                │  - Closest reference matching      │     │             │
                │  - Contamination detection         │     └─────────────┘
                └────────────────────────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐     ┌─────────────┐
                │  MODULE 4: Annotation              │────▶│ GFF/GBK     │
                │  - Gene prediction                 │     │ Annotations │
                │  - Functional annotation           │     │             │
                │  - Specialized databases           │     └─────────────┘
                └────────────────────────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────────────────────────┐
                │  MODULE 5: Specialized Analysis (Parallel Branches)    │
                │                                                        │
                │  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐  │
                │  │   AMR    │ │Virulence │ │  Typing  │ │ Plasmids │  │
                │  │ Analysis │ │ Factors  │ │MLST/cgMLST│ │ & MGEs   │  │
                │  └──────────┘ └──────────┘ └──────────┘ └──────────┘  │
                │                                                        │
                │  ┌──────────┐ ┌──────────┐ ┌──────────┐ ┌──────────┐  │
                │  │ Prophage │ │ CRISPR   │ │ Genomic  │ │ Secretion│  │
                │  │Detection │ │ Arrays   │ │ Islands  │ │ Systems  │  │
                │  └──────────┘ └──────────┘ └──────────┘ └──────────┘  │
                └────────────────────────────────────────────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐     ┌─────────────┐
                │  MODULE 6: Anomaly Detection       │────▶│ Anomaly     │
                │  - Foreign element identification  │     │ Report      │
                │  - Mosaic genome detection         │     │             │
                │  - Horizontal gene transfer        │     └─────────────┘
                │  - Synthetic biology signatures    │
                └────────────────────────────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐     ┌─────────────┐
                │  MODULE 7: Phylogenetic Analysis   │────▶│ Phylogeny   │
                │  - (Reserved for future use)       │     │ Trees       │
                └────────────────────────────────────┘     └─────────────┘
                                 │
                                 ▼
                ┌────────────────────────────────────┐     ┌─────────────┐
                │  MODULE 8: Report Generation       │────▶│ HTML Report │
                │  - Aggregate all results           │     │ with        │
                │  - Generate visualizations         │     │ Interactive │
                │  - Create clinical summary         │     │ Visuals     │
                │  - Export annotations              │     └─────────────┘
                └────────────────────────────────────┘
```

## Directory Structure

```
genome_pipeline/
│
├── Snakefile                    # Main entry point
├── config/
│   ├── config.yaml              # Main configuration
│   ├── databases.yaml           # Database paths
│   └── test_sample.yaml         # Test configuration
│
├── workflow/
│   ├── rules/
│   │   ├── common.smk           # Shared rules and functions
│   │   ├── qc.smk               # Module 1: Quality Control
│   │   ├── assembly.smk         # Module 2: Assembly
│   │   ├── identification.smk   # Module 3: Identification
│   │   ├── annotation.smk       # Module 4: Annotation
│   │   ├── amr.smk              # Module 5a: AMR Analysis
│   │   ├── virulence.smk        # Module 5b: Virulence Factors
│   │   ├── mobile_elements.smk  # Module 5d: Plasmids/MGEs
│   │   ├── prophages.smk        # Module 5e: Prophage Detection
│   │   ├── anomaly.smk          # Module 6: Anomaly Detection
│   │   └── report.smk           # Module 8: Report Generation
│   │
│   └── envs/                    # Conda environments
│       ├── qc.yaml
│       ├── assembly.yaml
│       ├── annotation.yaml
│       └── ...
│
├── scripts/                     # Helper scripts
│   ├── api_helpers/             # API-based analysis helpers
│   └── generate_report_direct.py
│
├── docs/                        # Documentation
│   ├── preparation_checklist.md
│   └── pipeline_plan.md
│
├── preparation_scripts/         # Setup scripts
│   ├── download_databases.sh
│   └── validate_installation.sh
│
└── output/                      # Results (created on run)
    └── {sample_id}/
        ├── 00_logs/
        ├── 01_qc/
        ├── 02_assembly/
        ├── 03_identification/
        ├── 04_annotation/
        ├── 05_specialized/
        ├── 06_anomaly/
        └── 08_report/
```

## Quick Start

```bash
# Basic usage - auto mode
snakemake --configfile config/config.yaml \
  --config sample_id=my_sample \
           reads.short.r1=reads_R1.fastq.gz \
           reads.short.r2=reads_R2.fastq.gz \
  --cores 8 --use-conda

# Long-read only
snakemake --configfile config/config.yaml \
  --config sample_id=my_sample \
           reads.long=reads_nanopore.fastq.gz \
  --cores 8 --use-conda

# Hybrid mode
snakemake --configfile config/config.yaml \
  --config sample_id=my_sample \
           reads.short.r1=reads_R1.fastq.gz \
           reads.short.r2=reads_R2.fastq.gz \
           reads.long=reads_nanopore.fastq.gz \
  --cores 8 --use-conda
```

---

*Document Version: 1.0*
*Created: December 2024*

