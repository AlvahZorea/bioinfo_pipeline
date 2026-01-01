# Genome Analysis Pipeline

A modular Snakemake pipeline for bacterial and viral genome analysis.

## Features

- **Quality Control**: FastQC, fastp, Porechop, NanoFilt, Kraken2 contamination detection
- **Assembly**: SPAdes (short), Flye (long), Unicycler (hybrid), with optional polishing
- **Identification**: Kraken2 taxonomy, fastANI, geNomad (virus/plasmid), CheckV
- **Annotation**: Bakta (bacteria), Pharokka (viruses/phages)
- **AMR Detection**: AMRFinderPlus
- **Virulence Factors**: BLAST vs VFDB
- **Mobile Elements**: MOB-suite
- **Prophage Detection**: PhiSpy
- **Report Generation**: Self-contained HTML with visualizations

## Requirements

- Snakemake >= 7.0
- Conda/Mamba
- Required databases (see `config/databases.yaml`)

## Quick Start

### 1. Configure Database Paths

Edit `config/databases.yaml` with paths to your databases:

```yaml
databases:
  kraken2:
    standard: "/path/to/kraken2/standard"
  bakta: "/path/to/bakta/db"
  # ... etc
```

### 2. Run the Pipeline

**With test data:**
```bash
snakemake --configfile config/test_sample.yaml --cores 8 --use-conda
```

**With custom sample:**
```bash
snakemake --configfile config/config.yaml \
    --config sample_id=my_sample \
             reads.short.r1=my_reads_R1.fastq.gz \
             reads.short.r2=my_reads_R2.fastq.gz \
    --cores 8 --use-conda
```

**Hybrid assembly (short + long reads):**
```bash
snakemake --configfile config/config.yaml \
    --config sample_id=my_sample \
             reads.short.r1=illumina_R1.fastq.gz \
             reads.short.r2=illumina_R2.fastq.gz \
             reads.long=nanopore.fastq.gz \
    --cores 8 --use-conda
```

### 3. Dry Run (Preview)

```bash
snakemake --configfile config/test_sample.yaml --cores 8 -n
```

## Directory Structure

```
genome_pipeline/
├── Snakefile                 # Main entry point
├── config/
│   ├── config.yaml           # Main configuration
│   ├── databases.yaml        # Database paths
│   └── test_sample.yaml      # Test configuration
├── workflow/
│   ├── rules/                # Snakemake rules
│   │   ├── common.smk        # Shared functions
│   │   ├── qc.smk            # Quality control
│   │   ├── assembly.smk      # Assembly
│   │   ├── identification.smk # Identification
│   │   ├── annotation.smk    # Annotation
│   │   ├── amr.smk           # AMR detection
│   │   ├── virulence.smk     # Virulence factors
│   │   ├── mobile_elements.smk # MGE detection
│   │   ├── prophages.smk     # Prophage detection
│   │   └── report.smk        # Report generation
│   └── envs/                 # Conda environments
└── output/                   # Results (created on run)
    └── {sample_id}/
        ├── 00_logs/
        ├── 01_qc/
        ├── 02_assembly/
        ├── 03_identification/
        ├── 04_annotation/
        ├── 05_specialized/
        └── 08_report/
```

## Output Files

| Module | Key Outputs |
|--------|-------------|
| QC | `01_qc/multiqc_report.html`, clean reads |
| Assembly | `02_assembly/contigs.fasta`, `quast/report.html` |
| Identification | `03_identification/taxonomy.json` |
| Annotation | `04_annotation/annotation.gff`, `proteins.faa` |
| AMR | `05_specialized/amr/amrfinder_results.tsv` |
| Virulence | `05_specialized/virulence/vfdb_results.tsv` |
| MGE | `05_specialized/mge/contig_report.txt` |
| Prophages | `05_specialized/prophages/prophage_coordinates.tsv` |
| Report | `08_report/genome_report.html` |

## Configuration Options

### Main Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `sample_id` | Sample identifier | "test_sample" |
| `mode` | Read mode: auto, short, long, hybrid | "auto" |
| `threads` | Number of threads | 8 |
| `modules.*` | Enable/disable modules | true |

### QC Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `qc.min_quality` | Minimum quality score | 20 |
| `qc.min_length` | Minimum read length | 50 |
| `qc.remove_contaminants` | Remove contaminated reads | true |
| `qc.contaminant_db` | Custom Kraken2 DB for contaminants | null |

### Assembly Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `assembly.assembler` | Assembler: auto, spades, flye, unicycler | "auto" |
| `assembly.min_contig_length` | Filter contigs below this length | 500 |
| `assembly.polish` | Enable polishing | true |

## Required Databases

See `preparation_scripts/download_databases.sh` for download instructions.

| Database | Size | Purpose |
|----------|------|---------|
| Kraken2 Standard | ~100 GB | Taxonomic classification |
| Bakta | ~30 GB | Bacterial annotation |
| geNomad | ~3 GB | Virus/plasmid detection |
| CheckV | ~2 GB | Viral quality assessment |
| AMRFinder | ~500 MB | AMR detection |
| VFDB | ~100 MB | Virulence factors |
| Pharokka | ~8 GB | Phage annotation |

## License

MIT License

