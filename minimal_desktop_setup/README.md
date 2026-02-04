# Minimal Desktop Testing Setup

This folder contains scripts for setting up a **lightweight testing environment** on your desktop/laptop.

## Contents

| File | Purpose |
|------|---------|
| `minimal_tools_install.sh` | Install essential bioinformatics tools via conda |
| `download_minimal_databases.sh` | Download ~15-25 GB of databases (vs ~350 GB full) |
| `verify_minimal_setup.sh` | Verify everything is working |
| `SERVER_VS_DESKTOP_PLANNING.md` | Worksheet to plan what runs where |

## Quick Start

```bash
# 1. Make scripts executable
chmod +x *.sh

# 2. Install tools (requires conda/mamba)
./minimal_tools_install.sh

# 3. Activate environment
conda activate genome-pipeline

# 4. Download minimal databases
./download_minimal_databases.sh ~/genome_pipeline_db

# 5. Verify setup
./verify_minimal_setup.sh ~/genome_pipeline_db
```

## What's Included vs Full Setup

| Component | Minimal (Desktop) | Full (Server) |
|-----------|-------------------|---------------|
| **Total Size** | ~15-25 GB | ~350 GB |
| **Kraken2 DB** | Mini (8 GB) | Standard (100 GB) |
| **Bakta DB** | Light (2 GB) | Full (30 GB) |
| **GTDB-Tk** | ❌ Not included | ✓ (85 GB) |
| **CheckV** | Optional (2 GB) | ✓ (2 GB) |
| **Pharokka** | ❌ Not included | ✓ (8 GB) |
| **InterProScan** | ❌ Not included | ✓ (50 GB) |

## What You Can Test

With the minimal setup you can:
- ✅ Assemble bacterial genomes (SPAdes, Flye)
- ✅ Basic species identification (Kraken2 mini)
- ✅ Gene annotation (Bakta light)
- ✅ AMR gene detection (AMRFinderPlus)
- ✅ MLST typing
- ✅ Virulence factor screening (VFDB)
- ✅ Quality control (FastQC, QUAST)

## What Requires Server

These analyses need the full database setup:
- ⚠️ Accurate taxonomy for rare organisms
- ⚠️ GTDB-Tk bacterial classification
- ⚠️ Phage annotation (Pharokka)
- ⚠️ Comprehensive protein domain analysis
- ⚠️ Processing many samples in parallel

## Next Steps

1. Fill out `SERVER_VS_DESKTOP_PLANNING.md` to inventory what you have
2. Decide on a workflow strategy
3. Start building the pipeline!
