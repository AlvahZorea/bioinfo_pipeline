# API Helper Scripts

These standalone scripts run bioinformatics analyses using web APIs instead of local databases.
**They are NOT part of the Snakemake pipeline** - use them for testing on a laptop without databases.

## Quick Start

```bash
# Run all analyses on a sample
python scripts/api_helpers/run_all_analyses.py output/citrobacter_sample --max-seqs 50

# Then regenerate the report
python scripts/generate_report_direct.py
```

## Individual Scripts

### 1. Extract Proteins
Extract protein sequences from Bakta GenBank annotation:
```bash
python scripts/api_helpers/extract_proteins.py bakta.gbff proteins.faa
```

### 2. NCBI BLAST for Per-Gene Taxonomy
Submit proteins to NCBI BLAST to classify each gene's taxonomic origin:
```bash
python scripts/api_helpers/ncbi_blast_api.py proteins.faa gene_taxonomy.json --max-seqs 100
```
⚠️ **Slow!** ~3 seconds per protein due to API rate limits.

### 3. ResFinder for AMR
For AMR gene detection, use the CGE web service:
```bash
# Show instructions for web submission
python scripts/api_helpers/resfinder_api.py contigs.fasta amr_summary.json

# Parse downloaded results
python scripts/api_helpers/resfinder_api.py contigs.fasta amr_summary.json --resfinder-output results.txt
```

### 4. VFDB for Virulence Factors
For virulence factor detection:
```bash
# Show instructions for web submission
python scripts/api_helpers/vfdb_api.py proteins.faa virulence_summary.json --web-only

# Or run NCBI BLAST with virulence filter (slower, less accurate)
python scripts/api_helpers/vfdb_api.py proteins.faa virulence_summary.json --max-seqs 50
```

## Web Services

| Analysis | Web Service | URL |
|----------|-------------|-----|
| AMR | ResFinder | https://cge.food.dtu.dk/services/ResFinder/ |
| Virulence | VFDB VFanalyzer | http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi |
| Plasmids | PlasmidFinder | https://cge.food.dtu.dk/services/PlasmidFinder/ |
| MLST | MLST | https://cge.food.dtu.dk/services/MLST/ |

## Output Format

All scripts output JSON files compatible with the pipeline format, which can be
placed in the sample output directory and used by `generate_report_direct.py`.

## Requirements

```bash
pip install biopython requests
```

