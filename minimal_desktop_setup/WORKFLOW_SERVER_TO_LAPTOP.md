# Workflow: Server → Laptop

This document describes the recommended workflow for your setup:
- **Server**: SPAdes, Flye, Kraken2, Prokka, GTDB-Tk, AMR tools
- **Laptop**: 8-16 GB RAM, ~20 GB disk

---

## Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│                           SERVER                                     │
│                                                                      │
│   Raw Reads (FASTQ)                                                  │
│        │                                                             │
│        ▼                                                             │
│   ┌─────────┐     ┌─────────┐     ┌─────────┐                       │
│   │ SPAdes  │────▶│ Kraken2 │────▶│ Prokka  │                       │
│   │Assembly │     │Taxonomy │     │Annotate │                       │
│   └─────────┘     └─────────┘     └─────────┘                       │
│        │               │               │                             │
│        ▼               ▼               ▼                             │
│   contigs.fasta   kraken.report   annotation.gff                    │
│                                                                      │
└────────────────────────┬────────────────────────────────────────────┘
                         │
                    Transfer via USB/network
                         │
                         ▼
┌─────────────────────────────────────────────────────────────────────┐
│                           LAPTOP                                     │
│                                                                      │
│   Import: contigs.fasta, kraken.report, annotation.gff              │
│        │                                                             │
│        ▼                                                             │
│   ┌──────────┐  ┌──────────┐  ┌──────────┐  ┌──────────┐           │
│   │AMRFinder │  │   MLST   │  │   VFDB   │  │  Report  │           │
│   │ (local)  │  │ (local)  │  │ (local)  │  │Generator │           │
│   └──────────┘  └──────────┘  └──────────┘  └──────────┘           │
│        │              │             │              │                 │
│        └──────────────┴─────────────┴──────────────┘                 │
│                              │                                       │
│                              ▼                                       │
│                    Final HTML Report                                 │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Step-by-Step Workflow

### PHASE 1: Run on Server

#### 1.1 Assembly (SPAdes)

```bash
# On SERVER
spades.py \
    -1 sample_R1.fastq.gz \
    -2 sample_R2.fastq.gz \
    -o assembly_output \
    --careful \
    -t 16

# Output needed: assembly_output/contigs.fasta
```

#### 1.2 Taxonomy (Kraken2)

```bash
# On SERVER
kraken2 \
    --db /path/to/kraken2_standard \
    --paired sample_R1.fastq.gz sample_R2.fastq.gz \
    --report kraken_report.txt \
    --output kraken_output.txt \
    --threads 16

# Output needed: kraken_report.txt
```

#### 1.3 Annotation (Prokka)

```bash
# On SERVER
prokka \
    --outdir prokka_output \
    --prefix sample \
    --cpus 16 \
    assembly_output/contigs.fasta

# Outputs needed: 
#   prokka_output/sample.gff
#   prokka_output/sample.faa (proteins)
#   prokka_output/sample.ffn (genes)
```

#### 1.4 (Optional) GTDB-Tk

```bash
# On SERVER (needs 64+ GB RAM)
gtdbtk classify_wf \
    --genome_dir genomes/ \
    --out_dir gtdbtk_output \
    --cpus 16

# Output needed: gtdbtk_output/classify/gtdbtk.bac120.summary.tsv
```

### PHASE 2: Transfer Files

**Files to transfer to laptop:**

```
sample_results/
├── assembly/
│   └── contigs.fasta          # ~5-50 MB per sample
├── taxonomy/
│   └── kraken_report.txt      # ~1 MB
├── annotation/
│   ├── sample.gff             # ~5-20 MB
│   ├── sample.faa             # ~2-10 MB
│   └── sample.ffn             # ~5-15 MB
└── gtdbtk/
    └── gtdbtk.summary.tsv     # ~1 KB
```

**Total per sample: ~20-100 MB** (very manageable!)

Transfer methods:
- USB drive
- `scp` / `rsync` over network
- Cloud sync (if allowed)

### PHASE 3: Run on Laptop

#### 3.1 AMR Detection

```bash
# On LAPTOP
conda activate gp-minimal

# AMRFinderPlus
amrfinder \
    --protein sample.faa \
    --nucleotide contigs.fasta \
    --gff sample.gff \
    --organism "Escherichia" \
    --output amr_results.tsv

# ABRicate (multiple databases)
abricate --db resfinder contigs.fasta > abricate_resfinder.tsv
abricate --db card contigs.fasta > abricate_card.tsv
abricate --db vfdb contigs.fasta > abricate_vfdb.tsv
```

#### 3.2 MLST Typing

```bash
# On LAPTOP
mlst contigs.fasta > mlst_results.tsv
```

#### 3.3 Virulence Factors

```bash
# On LAPTOP (using downloaded VFDB)
blastn \
    -query contigs.fasta \
    -db ~/genome_db/vfdb/vfdb_nucl \
    -outfmt "6 qseqid sseqid pident length qlen slen evalue" \
    -evalue 1e-10 \
    -out vfdb_results.tsv
```

#### 3.4 Generate Report

```bash
# On LAPTOP
# The pipeline will aggregate all results into HTML report
snakemake --configfile config.yaml report
```

---

## File Organization

### On Server

```
/data/projects/my_analysis/
├── raw_reads/
│   ├── sample1_R1.fastq.gz
│   └── sample1_R2.fastq.gz
├── assembly/
│   └── sample1/contigs.fasta
├── taxonomy/
│   └── sample1/kraken_report.txt
└── annotation/
    └── sample1/
        ├── sample1.gff
        ├── sample1.faa
        └── sample1.ffn
```

### On Laptop

```
~/genome_analysis/
├── genome_db/                    # Databases (~3-4 GB)
│   ├── amrfinder/
│   ├── vfdb/
│   └── config.yaml
├── pipeline/                     # Pipeline code
│   └── (Snakemake files)
├── samples/                      # Transferred from server
│   └── sample1/
│       ├── contigs.fasta
│       ├── kraken_report.txt
│       ├── sample1.gff
│       └── sample1.faa
└── results/                      # Laptop analysis outputs
    └── sample1/
        ├── amr_results.tsv
        ├── mlst_results.tsv
        ├── vfdb_results.tsv
        └── report.html
```

---

## Quick Reference: What Runs Where

| Analysis | Server | Laptop | Notes |
|----------|--------|--------|-------|
| **Assembly** | ✅ | ❌ | SPAdes needs RAM |
| **Kraken2** | ✅ | ❌ | Standard DB too large |
| **GTDB-Tk** | ✅ | ❌ | Needs 64+ GB RAM |
| **Prokka** | ✅ | ⚠️ | Can run on laptop but slower |
| **AMRFinderPlus** | ✅ | ✅ | Small DB, fast |
| **MLST** | ✅ | ✅ | Small DB, instant |
| **ABRicate** | ✅ | ✅ | Built-in DBs |
| **QUAST** | ✅ | ✅ | No DB needed |
| **FastQC** | ✅ | ✅ | No DB needed |
| **Report Gen** | ⚠️ | ✅ | Best on laptop for iteration |

---

## Example: Complete Sample Analysis

```bash
#=============================================================================
# SERVER COMMANDS (run first)
#=============================================================================

# Variables
SAMPLE="sample1"
READS_DIR="/data/reads"
OUT_DIR="/data/analysis/${SAMPLE}"

mkdir -p ${OUT_DIR}/{assembly,taxonomy,annotation}

# 1. Assembly
spades.py \
    -1 ${READS_DIR}/${SAMPLE}_R1.fastq.gz \
    -2 ${READS_DIR}/${SAMPLE}_R2.fastq.gz \
    -o ${OUT_DIR}/assembly \
    --careful -t 16

# 2. Kraken2
kraken2 \
    --db /databases/kraken2_standard \
    --paired ${READS_DIR}/${SAMPLE}_R1.fastq.gz ${READS_DIR}/${SAMPLE}_R2.fastq.gz \
    --report ${OUT_DIR}/taxonomy/kraken_report.txt \
    --output ${OUT_DIR}/taxonomy/kraken_output.txt

# 3. Prokka
prokka \
    --outdir ${OUT_DIR}/annotation \
    --prefix ${SAMPLE} \
    --cpus 16 \
    ${OUT_DIR}/assembly/contigs.fasta

# 4. Package for transfer
tar -czf ${SAMPLE}_for_laptop.tar.gz \
    ${OUT_DIR}/assembly/contigs.fasta \
    ${OUT_DIR}/taxonomy/kraken_report.txt \
    ${OUT_DIR}/annotation/${SAMPLE}.gff \
    ${OUT_DIR}/annotation/${SAMPLE}.faa \
    ${OUT_DIR}/annotation/${SAMPLE}.ffn

#=============================================================================
# LAPTOP COMMANDS (after transfer)
#=============================================================================

# Variables
SAMPLE="sample1"
SAMPLE_DIR=~/genome_analysis/samples/${SAMPLE}

# Extract transferred files
mkdir -p ${SAMPLE_DIR}
tar -xzf ${SAMPLE}_for_laptop.tar.gz -C ${SAMPLE_DIR}

# Activate environment
conda activate gp-minimal

# 5. AMR Detection
amrfinder \
    -p ${SAMPLE_DIR}/${SAMPLE}.faa \
    -n ${SAMPLE_DIR}/contigs.fasta \
    -g ${SAMPLE_DIR}/${SAMPLE}.gff \
    -o ${SAMPLE_DIR}/amr_results.tsv

# 6. MLST
mlst ${SAMPLE_DIR}/contigs.fasta > ${SAMPLE_DIR}/mlst_results.tsv

# 7. Virulence
abricate --db vfdb ${SAMPLE_DIR}/contigs.fasta > ${SAMPLE_DIR}/vfdb_results.tsv

# 8. Generate report (once pipeline is built)
# snakemake --config sample=${SAMPLE} report
```

---

## Troubleshooting

### "Not enough disk space on laptop"
- Delete intermediate files after analysis
- Process one sample at a time
- Keep only final reports, archive rest to external drive

### "Analysis too slow on laptop"
- Reduce thread count if RAM is limited
- Consider running on server instead
- Use `--fast` modes where available

### "Different results between server and laptop AMR"
- Ensure same database versions
- Check AMRFinderPlus database date
- Server results are canonical; laptop is supplementary
