# Bioinformatics Pipeline - Project Summary

**Last Updated:** 2024  
**Status:** Active Development - Core Pipeline Complete, Testing Phase  
**Repository:** https://github.com/AlvahZorea/bioinfo_pipeline

---

## Table of Contents

1. [Project Overview](#project-overview)
2. [Architecture & Design Philosophy](#architecture--design-philosophy)
3. [Development History](#development-history)
4. [Pipeline Modules](#pipeline-modules)
5. [Current State Assessment](#current-state-assessment)
6. [Known Limitations & Challenges](#known-limitations--challenges)
7. [Next Steps](#next-steps)
8. [Usage Guide](#usage-guide)
9. [Technical Details for LLMs](#technical-details-for-llms)

---

## Project Overview

### Purpose

This project implements a comprehensive **Snakemake-based bioinformatics pipeline** for bacterial and viral genome analysis with a specific focus on **detecting genomic editing events and anomalies**. The pipeline is designed to identify:

- Foreign gene insertions (horizontal gene transfer, synthetic biology)
- Structural variants (deletions, insertions, rearrangements)
- Integrated mobile elements (plasmids, ICEs)
- Chimeric/mosaic genomes
- Novel or synthetic genes
- Regulatory element insertions

### Key Requirements

1. **Modular Design**: Each analysis module can be enabled/disabled independently
2. **Flexible Input**: Supports short-read (Illumina), long-read (Nanopore), and hybrid assemblies
3. **Comprehensive Analysis**: QC → Assembly → Identification → Annotation → Specialized Analysis → Anomaly Detection → Reporting
4. **Anomaly Detection**: Core focus on detecting genomic modifications and foreign elements
5. **Self-Contained Reports**: HTML reports with interactive visualizations

### Target Use Cases

- **Quality Control**: Standard genome assembly and annotation workflows
- **Pathogen Analysis**: AMR, virulence factor detection
- **Forensic Genomics**: Detection of engineered/modified genomes
- **Research**: Identification of horizontal gene transfer events
- **Surveillance**: Monitoring for novel genetic elements

---

## Architecture & Design Philosophy

### Technology Stack

- **Workflow Engine**: Snakemake (v7+)
- **Environment Management**: Conda/Mamba
- **Language**: Python 3.8+ (with BioPython, Plotly)
- **Visualization**: Plotly.js (embedded in HTML reports)
- **Version Control**: Git/GitHub

### Design Principles

1. **Modularity**: Each module is a separate `.smk` file with isolated dependencies
2. **Configurability**: All parameters in `config/config.yaml`, database paths in `config/databases.yaml`
3. **Reproducibility**: Conda environments for each module ensure consistent tool versions
4. **Flexibility**: Can run full pipeline or individual modules
5. **Extensibility**: Easy to add new analysis modules or tools

### Directory Structure

```
bioinfo_pipeline/
├── config/
│   ├── config.yaml          # Main pipeline configuration
│   └── databases.yaml        # Database paths
├── workflow/
│   ├── rules/                # Snakemake rule files (.smk)
│   │   ├── qc.smk
│   │   ├── assembly.smk
│   │   ├── identification.smk
│   │   ├── annotation.smk
│   │   ├── amr.smk
│   │   ├── virulence.smk
│   │   ├── mobile_elements.smk
│   │   ├── prophages.smk
│   │   ├── anomaly.smk       # Core anomaly detection
│   │   ├── report.smk
│   │   └── common.smk        # Helper functions
│   └── envs/                  # Conda environment files
│       ├── qc.yaml
│       ├── assembly.yaml
│       └── ...
├── scripts/                  # Helper Python scripts
│   ├── api_helpers/          # API-based analysis scripts
│   └── ...
├── Snakefile                 # Main entry point
└── output/                   # Pipeline outputs (per sample)
    └── {sample_id}/
        ├── 01_qc/
        ├── 02_assembly/
        ├── 03_identification/
        ├── 04_annotation/
        ├── 05_specialized/
        ├── 06_anomaly/
        └── 08_report/
```

---

## Development History

### Phase 1: Initial Planning & Design (Pre-Implementation)

**Goal**: Establish pipeline architecture and tool selection

**Key Decisions**:
- Chose Snakemake over Nextflow (better Python integration, simpler syntax)
- Decided on modular structure with separate rule files per module
- Selected tools based on user requirements and expert recommendations:
  - **QC**: FastQC, fastp, Kraken2 (contamination), MultiQC
  - **Assembly**: SPAdes (short), Flye (long), Unicycler (hybrid)
  - **Identification**: Kraken2, fastANI, geNomad, CheckV
  - **Annotation**: Bakta (bacterial), Pharokka (viral)
  - **AMR**: AMRFinderPlus
  - **Virulence**: VFDB (BLAST)
  - **MGEs**: MOB-suite
  - **Prophages**: PhiSpy

**Outcome**: Created `genome_pipeline_plan.md.txt` with comprehensive design document

### Phase 2: Core Pipeline Implementation

**Goal**: Build basic pipeline structure and essential modules

**Accomplishments**:
- Created Snakemake directory structure
- Implemented all 8 core modules (QC through Report)
- Set up Conda environment definitions
- Created configuration system (`config.yaml`, `databases.yaml`)
- Implemented helper functions in `common.smk`
- Built basic HTML report generation

**Challenges**:
- Handling optional inputs (not all modules always run)
- Configuring database paths flexibly
- Making report generation work with partial data

**Outcome**: Functional pipeline skeleton ready for testing

### Phase 3: Mock Data & Report Testing

**Goal**: Test report generation without running full pipeline

**Accomplishments**:
- Created mock data structure (`mock_data/ecoli_salmonella_insert/`)
- Implemented `target_report_only` mode for report-only runs
- Generated test HTML reports with visualizations
- Tested with real `citrobacter` sample data from server outputs

**Key Insight**: Need to support "report-only" mode for using existing analysis outputs

**Outcome**: Report generation validated, visualization working

### Phase 4: Anomaly Detection Module

**Goal**: Implement core anomaly detection capabilities

**Initial Implementation**:
- Mash-based reference selection
- MUMmer (nucmer) alignment for structural variants
- Per-gene BLAST taxonomy analysis
- Basic foreign gene detection

**Visualization**:
- Circular genome maps with foreign gene highlighting
- Linear alignment views
- Anomaly alert sections in reports

**Outcome**: Basic anomaly detection functional

### Phase 5: API Helper Scripts (Laptop Testing)

**Goal**: Enable analysis on laptop without large databases

**Accomplishments**:
- Created `scripts/api_helpers/` for web API-based analysis:
  - NCBI BLAST API for per-gene taxonomy
  - ResFinder API instructions
  - VFDB API instructions
- Handled SSL certificate issues on Windows
- Created master orchestration script

**Limitation**: API calls are slower and have rate limits; not ideal for production

**Outcome**: Enabled testing on laptop, but identified need for server access

### Phase 6: Table Exercise & Blind Spot Analysis

**Goal**: Identify pipeline limitations through theoretical scenarios

**Scenarios Analyzed**:
1. Hybrid/mosaic genome (half E. coli, half Salmonella)
2. E. coli with integrated plasmid
3. E. coli with free Salmonella plasmid
4. E. coli with integrated Salmonella plasmid
5. Gene location swaps
6. Gene deletions
7. Unknown origin gene insertions
8. Promoter insertions

**Identified Blind Spots**:
- **Synteny analysis**: No gene order comparison
- **Intergenic regions**: Not scanned for foreign elements
- **Integrated elements**: No detection of plasmid backbone genes in chromosome
- **Deletion annotation**: Gaps detected but not mapped to genes
- **Novel genes**: No flagging of genes with no BLAST hits
- **Minority species**: No alert for mixed taxonomy
- **Foreign gene clustering**: Adjacent foreign genes not grouped

**Outcome**: Comprehensive list of improvements needed

### Phase 7: Blind Spot Fixes (Current)

**Goal**: Implement all identified improvements

**Implementations**:

1. **NCBI Reference Auto-Download**
   - Rule: `download_reference_ncbi`
   - Uses NCBI Datasets API to fetch reference based on Kraken2 taxonomy
   - Caches downloads for reuse

2. **Synteny Analysis**
   - Rule: `synteny_analysis`
   - Compares gene order between assembly and reference GFFs
   - Detects rearrangements, inversions, translocations

3. **Intergenic Region Scanning**
   - Rules: `extract_intergenic`, `blast_intergenic`
   - Extracts regions >100bp between genes
   - BLASTs against nt database to find foreign regulatory elements

4. **Integrated Element Detection**
   - Rule: `detect_integrated_elements`
   - Scans for plasmid backbone genes (repA, traA, mobA, integrase) in chromosome
   - Clusters adjacent hits to identify integrated plasmids/ICEs

5. **Deletion Annotation**
   - Enhanced `nucmer_align` rule
   - Maps alignment gaps to reference GFF
   - Reports which genes are deleted

6. **Novel Gene Flagging**
   - Enhanced `blast_genes` rule
   - Flags genes with no hits or <50% identity as potentially synthetic/novel

7. **Minority Species Alert**
   - Enhanced `parse_kraken_taxonomy` rule
   - Flags non-dominant genera >0.5% threshold

8. **Foreign Gene Clustering**
   - Enhanced `aggregate_anomalies` rule
   - Clusters adjacent foreign genes (within 10kb)
   - Flags 3+ gene clusters as integrated elements

9. **Report Enhancements**
   - Severity-based anomaly alerts (low/medium/high)
   - Minority species warnings
   - Foreign gene cluster tables
   - Novel gene tables
   - Deletion gene listings
   - Synteny visualization
   - Intergenic element tables
   - Integrated element detection display

**Outcome**: Comprehensive anomaly detection system implemented

---

## Pipeline Modules

### Module 1: Quality Control (`qc.smk`)

**Purpose**: Filter and assess read quality, remove contaminants

**Tools**:
- **FastQC**: Read quality assessment
- **fastp**: Adapter trimming, quality filtering (short reads)
- **Porechop**: Adapter trimming (long reads)
- **NanoFilt/Filtlong**: Quality filtering (long reads)
- **Kraken2**: Contamination screening
- **MultiQC**: Aggregate QC reports

**Outputs**:
- `01_qc/multiqc_report.html`
- Filtered reads (if `remove_contaminants: true`)

**Key Features**:
- Modular contaminant removal (can specify known contaminant DB)
- Standard quality thresholds (configurable)
- Separate handling for short vs long reads

### Module 2: Assembly (`assembly.smk`)

**Purpose**: Assemble reads into contigs

**Tools**:
- **SPAdes**: Short-read assembly
- **Flye**: Long-read assembly
- **Unicycler**: Hybrid assembly
- **Medaka**: Nanopore polishing
- **Pilon**: Short-read polishing
- **QUAST**: Assembly quality assessment

**Outputs**:
- `02_assembly/contigs.fasta`
- `02_assembly/assembly_stats.json`
- `02_assembly/quast/report.html`

**Key Features**:
- Auto-selects assembler based on read type
- Optional polishing (default: enabled)
- Configurable minimum contig length
- Assembly statistics tracking

### Module 3: Identification (`identification.smk`)

**Purpose**: Taxonomic identification and organism type determination

**Tools**:
- **Kraken2**: Taxonomic classification
- **fastANI**: Average Nucleotide Identity to references
- **geNomad**: Virus/plasmid identification
- **CheckV**: Viral genome quality

**Outputs**:
- `03_identification/taxonomy.json`
- `03_identification/organism_type.txt`
- `03_identification/kraken2_report.txt`

**Key Features**:
- Minority species detection (>0.5% threshold)
- Automatic organism type determination (bacteria/virus/archaea)
- Top taxonomic hits tracking

### Module 4: Annotation (`annotation.smk`)

**Purpose**: Functional annotation of genes

**Tools**:
- **Bakta**: Bacterial genome annotation
- **Pharokka**: Viral genome annotation

**Outputs**:
- `04_annotation/bakta/{sample_id}.gff3`
- `04_annotation/bakta/{sample_id}.gbff`
- `04_annotation/bakta/{sample_id}.faa`
- `04_annotation/annotation_stats.json`

**Key Features**:
- Organism-specific annotation (Bakta for bacteria)
- Comprehensive functional annotation
- GFF3 output for downstream analysis

### Module 5: Specialized Analysis

#### 5a: AMR (`amr.smk`)

**Tool**: AMRFinderPlus  
**Outputs**: `05_specialized/amr/amr_summary.json`

#### 5b: Virulence (`virulence.smk`)

**Tool**: VFDB (BLAST)  
**Outputs**: `05_specialized/virulence/virulence_summary.json`

#### 5c: Mobile Elements (`mobile_elements.smk`)

**Tool**: MOB-suite  
**Outputs**: `05_specialized/mge/mge_summary.json`

#### 5d: Prophages (`prophages.smk`)

**Tool**: PhiSpy  
**Outputs**: `05_specialized/prophages/prophage_summary.json`

### Module 6: Anomaly Detection (`anomaly.smk`) ⭐ **CORE MODULE**

**Purpose**: Detect genomic editing events and foreign elements

**Analysis Components**:

1. **Reference Selection**
   - `mash_sketch`: Create Mash sketch of assembly
   - `select_reference`: Find closest reference using Mash
   - `download_reference_ncbi`: Auto-download reference from NCBI

2. **Reference Alignment**
   - `nucmer_align`: Align assembly to reference (MUMmer)
   - Detects insertions, deletions
   - Maps deletions to reference genes

3. **Per-Gene Taxonomy**
   - `extract_proteins`: Extract proteins from Bakta annotation
   - `blast_genes`: BLAST each gene against nr database
   - Identifies foreign genes (different genus than host)
   - Flags novel genes (no hits or <50% identity)

4. **Synteny Analysis**
   - `synteny_analysis`: Compare gene order (assembly vs reference)
   - Detects rearrangements, inversions, translocations

5. **Intergenic Scanning**
   - `extract_intergenic`: Extract intergenic regions >100bp
   - `blast_intergenic`: BLAST regions to find foreign regulatory elements

6. **Integrated Element Detection**
   - `detect_integrated_elements`: Scan for plasmid backbone genes in chromosome
   - Clusters adjacent hits to identify integrated plasmids/ICEs

7. **Aggregation**
   - `aggregate_anomalies`: Combines all results
   - Clusters adjacent foreign genes
   - Calculates anomaly severity (low/medium/high)

**Outputs**:
- `06_anomaly/anomaly_summary.json`
- `06_anomaly/reference/reference_info.json`
- `06_anomaly/alignment/alignment.json`
- `06_anomaly/gene_taxonomy/gene_taxonomy.json`
- `06_anomaly/synteny/synteny_report.json`
- `06_anomaly/intergenic/intergenic_taxonomy.json`
- `06_anomaly/integrated_elements/integrated_elements.json`

### Module 8: Report Generation (`report.smk`)

**Purpose**: Generate comprehensive HTML report

**Features**:
- Self-contained HTML (all assets embedded)
- Interactive Plotly visualizations:
  - Circular genome map (GC content, foreign genes)
  - Linear contig view with anomalies
  - Reference alignment visualization
  - Contig length distribution
  - Taxonomy pie chart
- Sections:
  - Executive summary
  - QC results
  - Assembly statistics
  - Taxonomic identification (with minority species alerts)
  - Annotation summary
  - **Anomaly detection details** (comprehensive)
  - AMR analysis
  - Virulence factors
  - Mobile elements
  - Prophages

**Outputs**:
- `08_report/genome_report.html`
- `08_report/all_results.json`

---

## Current State Assessment

### ✅ What's Working

1. **Pipeline Structure**: Complete modular architecture
2. **Core Modules**: All 8 modules implemented and functional
3. **Anomaly Detection**: Comprehensive system with 7 analysis components
4. **Report Generation**: Rich HTML reports with visualizations
5. **Configuration System**: Flexible config files
6. **Conda Integration**: Environment management working
7. **Mock Data Testing**: Report generation tested with artificial data
8. **Real Data Testing**: Report generated from `citrobacter` server outputs

### ⚠️ Known Limitations

1. **Database Dependencies**: Many tools require large databases not yet downloaded:
   - Kraken2 database
   - Bakta database
   - NCBI BLAST databases (nr, nt)
   - RefSeq Mash sketches
   - AMRFinderPlus database
   - VFDB database
   - MOB-suite databases
   - CheckV database
   - geNomad database

2. **NCBI Datasets API**: Reference auto-download requires `datasets` CLI tool installation

3. **API-Based Analysis**: Current API helpers are workarounds; not production-ready:
   - Rate limits
   - SSL certificate issues on Windows
   - Slower than local databases

4. **Testing Limitations**: 
   - Only tested with pre-existing outputs (not full pipeline run)
   - No validation on actual genome editing scenarios yet

5. **Missing Validations**:
   - No automated tests
   - No validation against known edited genomes
   - No performance benchmarking

### 📊 Code Quality

- **Structure**: Well-organized, modular
- **Documentation**: Inline comments, docstrings
- **Error Handling**: Basic (needs improvement)
- **Logging**: Basic (needs structured logging)

---

## Known Limitations & Challenges

### Technical Challenges

1. **Database Size**: Bioinformatics databases are massive (100s of GB)
   - **Solution**: Download on server, use symlinks or network storage

2. **Tool Installation**: Many tools require complex dependencies
   - **Solution**: Conda environments handle this, but need to be tested

3. **Reference Genome Selection**: Auto-download may not always find best reference
   - **Solution**: Allow user override in config

4. **BLAST Performance**: Per-gene BLAST is slow for large genomes
   - **Solution**: Parallelization, consider DIAMOND for faster searches

5. **Memory Requirements**: Some tools (SPAdes, Bakta) need significant RAM
   - **Solution**: Configurable memory limits, use appropriate assemblers

### Biological Challenges

1. **False Positives**: Foreign genes may be legitimate HGT
   - **Solution**: Context matters - clusters more suspicious than isolated genes

2. **Novel Genes**: No BLAST hit doesn't guarantee synthetic origin
   - **Solution**: Flag as "suspicious" not "definitely synthetic"

3. **Reference Quality**: Poor reference choice leads to poor alignment
   - **Solution**: Use Mash distance threshold, allow user override

4. **Mixed Samples**: Contamination can trigger false anomaly alerts
   - **Solution**: Minority species alerts help identify this

---

## Next Steps

### Immediate (Before Server Access)

1. **Documentation**
   - ✅ Create this project summary (DONE)
   - Create user guide with examples
   - Document database download procedures
   - Create troubleshooting guide

2. **Code Cleanup**
   - Review and fix any linting issues
   - Add error handling improvements
   - Add progress logging
   - Validate all file paths exist before use

3. **Configuration Validation**
   - Add config validation script
   - Check database paths exist
   - Verify tool installations

4. **Test Data Preparation**
   - Create test datasets for each scenario from table exercise
   - Prepare validation genomes (known edited vs unedited)

### Short-Term (With Server Access)

1. **Database Setup**
   - Download all required databases:
     ```
     - Kraken2 standard database (~100GB)
     - Bakta database (~10GB)
     - NCBI BLAST databases (nr, nt) (~200GB)
     - RefSeq Mash sketches (~50GB)
     - AMRFinderPlus database (~500MB)
     - VFDB database (~1GB)
     - MOB-suite databases (~5GB)
     - CheckV database (~10GB)
     - geNomad database (~20GB)
     ```
   - Update `config/databases.yaml` with paths
   - Test database accessibility

2. **Full Pipeline Testing**
   - Run complete pipeline on test samples
   - Test each module independently
   - Validate outputs match expected formats
   - Test error handling (missing files, tool failures)

3. **Performance Optimization**
   - Profile slow steps (BLAST, alignment)
   - Optimize parallelization
   - Consider DIAMOND for protein searches
   - Cache intermediate results where possible

4. **Validation Testing**
   - Test on known edited genomes (if available)
   - Test on unedited control genomes
   - Validate anomaly detection accuracy
   - Test edge cases (very small genomes, very large genomes)

### Medium-Term (Production Readiness)

1. **Automated Testing**
   - Create test suite (pytest)
   - Test each module with known inputs/outputs
   - Integration tests for full pipeline
   - Continuous integration (GitHub Actions)

2. **Error Handling & Logging**
   - Structured logging (JSON format)
   - Better error messages
   - Graceful degradation (continue with partial results)
   - Retry logic for transient failures

3. **Performance Monitoring**
   - Track execution time per module
   - Memory usage monitoring
   - Resource usage reports

4. **User Experience**
   - Better progress indicators
   - Clearer error messages
   - Usage examples
   - FAQ document

5. **Advanced Features**
   - Batch processing (multiple samples)
   - Resume from checkpoint
   - Dry-run mode improvements
   - Configuration templates for common use cases

### Long-Term (Enhancements)

1. **Additional Anomaly Detection Methods**
   - GC content/codon usage analysis (if needed)
   - K-mer based anomaly detection
   - Machine learning for anomaly classification

2. **Enhanced Visualizations**
   - Interactive genome browser
   - Comparison views (multiple samples)
   - Export to standard formats (PDF, PNG)

3. **Integration**
   - Galaxy workflow export
   - Nextflow conversion option
   - API for programmatic access

4. **Scalability**
   - Cluster execution (SLURM, SGE)
   - Cloud deployment options
   - Containerization (Docker, Singularity)

---

## Usage Guide

### Basic Usage

```bash
# Install Snakemake
pip install snakemake

# Run full pipeline
snakemake --configfile config/config.yaml --cores 8 --use-conda

# Run specific module only
snakemake --configfile config/config.yaml --cores 8 --use-conda 08_report/genome_report.html

# Dry run (see what would be executed)
snakemake --configfile config/config.yaml --cores 8 --use-conda -n

# Run with custom sample
snakemake --configfile config/config.yaml \
    --config sample_id=my_sample \
             reads.short.r1=reads_R1.fastq.gz \
             reads.short.r2=reads_R2.fastq.gz \
    --cores 8 --use-conda
```

### Configuration

Edit `config/config.yaml`:
- Set `sample_id` and `sample_name`
- Configure read paths (`reads.short.r1`, `reads.short.r2`, `reads.long`)
- Enable/disable modules (`modules.qc`, `modules.assembly`, etc.)
- Adjust parameters per module

Edit `config/databases.yaml`:
- Set paths to all required databases
- Use absolute paths or paths relative to project root

### Module-Specific Usage

**Report-only mode** (using existing outputs):
```yaml
# In config.yaml
target_report_only: true
output_dir: "path/to/existing/outputs"
sample_id: "existing_sample"
```

**Anomaly detection only**:
```bash
snakemake --configfile config/config.yaml \
    --config modules.qc=false modules.assembly=false \
             modules.identification=false modules.annotation=false \
             modules.amr=false modules.virulence=false \
             modules.mobile_elements=false modules.prophages=false \
    --cores 8 --use-conda 06_anomaly/anomaly_summary.json
```

### Output Structure

```
output/{sample_id}/
├── 00_logs/              # Pipeline logs
├── 01_qc/                # Quality control reports
├── 02_assembly/          # Assembly files and stats
├── 03_identification/     # Taxonomic identification
├── 04_annotation/        # Gene annotations (GFF, GBFF, FAA)
├── 05_specialized/       # AMR, virulence, MGEs, prophages
│   ├── amr/
│   ├── virulence/
│   ├── mge/
│   └── prophages/
├── 06_anomaly/           # Anomaly detection results
│   ├── reference/        # Reference selection
│   ├── alignment/       # MUMmer alignment
│   ├── gene_taxonomy/    # Per-gene BLAST results
│   ├── synteny/         # Gene order analysis
│   ├── intergenic/      # Intergenic region analysis
│   └── integrated_elements/  # Integrated element detection
└── 08_report/           # Final HTML report
    └── genome_report.html
```

---

## Technical Details for LLMs

### For Another Cursor Instance or LLM

If you're picking up this project, here's what you need to know:

#### Key Files to Understand

1. **`Snakefile`**: Main entry point
   - Defines global variables (`SAMPLE_ID`, `OUTPUT_DIR`, `THREADS`)
   - Includes all module rule files
   - Defines `all` rule (final targets)

2. **`workflow/rules/common.smk`**: Helper functions
   - `get_sample_id(config)`: Extract sample ID
   - `get_sample_output(config)`: Get output directory
   - `get_threads(config)`: Get thread count
   - `detect_read_mode(config)`: Auto-detect short/long/hybrid
   - `is_module_enabled(config, module)`: Check if module enabled
   - `get_database_path(config, db_name, default)`: Get DB path

3. **`config/config.yaml`**: Main configuration
   - Sample info, read paths, module flags
   - Per-module parameters
   - Anomaly detection thresholds

4. **`workflow/rules/anomaly.smk`**: Core anomaly detection
   - Most complex module
   - Multiple analysis components
   - Aggregates results into `anomaly_summary.json`

5. **`workflow/rules/report.smk`**: Report generation
   - Python functions generate HTML
   - Uses Plotly for visualizations
   - Aggregates all module outputs

#### Common Patterns

**Rule Structure**:
```python
rule rule_name:
    input: "path/to/input"
    output: "path/to/output"
    params: param_value
    threads: 8
    log: "path/to/log"
    conda: "../envs/module.yaml"
    shell: "command {input} {output}"
    # OR
    run:
        # Python code
```

**Optional Inputs**:
- Check if file exists: `if os.path.exists(path):`
- Use in `aggregate_results`: Check existence before loading

**Module Dependencies**:
- Anomaly module depends on: annotation (for GFF), identification (for taxonomy)
- Report depends on: all upstream modules (but handles missing gracefully)

#### Database Paths

Databases are configured in `config/databases.yaml`:
- Access via: `config.get("databases", {}).get("db_name", "")`
- Check existence before using: `if db_path and os.path.exists(db_path):`

#### Error Handling Pattern

```python
try:
    # Run tool
    result = subprocess.run(cmd, shell=True)
    if result.returncode != 0:
        # Handle error
except Exception as e:
    # Create empty output, log error
    open(output.file, 'w').close()
```

#### Testing Approach

1. **Mock Data**: Use `mock_data/` for testing report generation
2. **Real Data**: Use existing server outputs (parse with `scripts/parse_existing_outputs.py`)
3. **Full Pipeline**: Requires all databases and tools installed

#### Adding New Features

1. **New Module**: 
   - Create `workflow/rules/new_module.smk`
   - Create `workflow/envs/new_module.yaml`
   - Add to `Snakefile` includes
   - Add to `get_final_outputs()` in `Snakefile`

2. **New Anomaly Detection Method**:
   - Add rule to `anomaly.smk`
   - Add output to `aggregate_anomalies` input
   - Parse results in `aggregate_anomalies`
   - Add visualization to `report.smk`

3. **New Report Section**:
   - Add function to `report.smk` (e.g., `generate_new_section()`)
   - Call in `generate_html_report()`
   - Add data loading in `generate_report` rule

#### Current State Summary

- **Status**: Core pipeline complete, anomaly detection comprehensive
- **Testing**: Limited (mock data + one real sample)
- **Databases**: Not yet downloaded (need server access)
- **Next Priority**: Full pipeline testing with real data and databases

#### Common Issues & Solutions

1. **"No reads provided"**: Set `target_report_only: true` in config for report-only mode
2. **"Database not found"**: Update `config/databases.yaml` with correct paths
3. **"Tool not found"**: Install via Conda (`--use-conda` flag)
4. **SSL errors (Windows)**: See `scripts/api_helpers/ncbi_blast_api.py` for workaround

---

## Conclusion

This pipeline represents a comprehensive solution for bacterial/viral genome analysis with a strong focus on anomaly detection. The architecture is modular and extensible, making it easy to add new analysis methods or tools.

**Current Status**: Ready for full testing once databases are available on server.

**Next Milestone**: Complete end-to-end testing with real data and validation against known edited genomes.

**Long-term Vision**: Production-ready pipeline for forensic genomics and research applications.

---

*For questions or contributions, see the GitHub repository: https://github.com/AlvahZorea/bioinfo_pipeline*

