#!/bin/bash
#===============================================================================
# ULTRA-MINIMAL LAPTOP SETUP (~10-12 GB total)
#===============================================================================
# For laptops with limited disk space (~20 GB free)
# Assumes heavy analysis (assembly, Kraken2, GTDB-Tk) runs on SERVER
#
# This setup focuses on:
#   - Pipeline development and testing
#   - AMR/virulence analysis (small databases)
#   - MLST typing
#   - Report generation
#   - Integrating pre-computed server results
#===============================================================================

set -e

DB_ROOT="${1:-$HOME/genome_db}"

echo "==============================================================================="
echo " ULTRA-MINIMAL LAPTOP SETUP"
echo " Target: ~10-12 GB total (tools + databases)"
echo " Database location: ${DB_ROOT}"
echo "==============================================================================="
echo ""

# Check available space
AVAILABLE=$(df -BG ~ | tail -1 | awk '{print $4}' | tr -d 'G')
echo "Available disk space: ${AVAILABLE} GB"
if [ "$AVAILABLE" -lt 15 ]; then
    echo "WARNING: Less than 15 GB free. Proceed with caution."
fi
echo ""

#-------------------------------------------------------------------------------
# STEP 1: Create minimal conda environment (~5-6 GB)
#-------------------------------------------------------------------------------
echo "=== STEP 1: Installing Tools via Conda ==="
echo ""

if ! command -v conda &> /dev/null && ! command -v mamba &> /dev/null; then
    echo "ERROR: conda/mamba not found. Install miniforge first:"
    echo "  https://github.com/conda-forge/miniforge"
    exit 1
fi

CONDA_CMD=$(command -v mamba || command -v conda)

echo "Creating ultra-minimal environment..."
$CONDA_CMD create -n gp-minimal -y \
    -c conda-forge -c bioconda \
    python=3.11 \
    snakemake-minimal \
    biopython \
    pandas \
    numpy \
    plotly \
    jinja2 \
    pyyaml \
    ncbi-amrfinderplus \
    mlst \
    abricate \
    blast \
    fastqc \
    multiqc \
    quast

echo ""
echo "Tools installed. Activate with: conda activate gp-minimal"
echo ""

#-------------------------------------------------------------------------------
# STEP 2: Download minimal databases (~3-4 GB)
#-------------------------------------------------------------------------------
echo "=== STEP 2: Downloading Minimal Databases ==="
echo ""

mkdir -p "${DB_ROOT}"/{amrfinder,mlst,vfdb,abricate}

# Activate environment for database downloads
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate gp-minimal

#--- AMRFinderPlus (~500 MB) ---
echo "Downloading AMRFinderPlus database..."
if [ ! -d "${DB_ROOT}/amrfinder/latest" ]; then
    amrfinder_update -d "${DB_ROOT}/amrfinder" || echo "AMRFinder DB failed"
else
    echo "  Already exists, skipping."
fi

#--- MLST (~500 MB - 1 GB) ---
echo "Updating MLST schemes..."
mlst --check 2>/dev/null || mlst --update 2>/dev/null || echo "MLST update failed"

#--- ABRicate databases (~200 MB total) ---
echo "Setting up ABRicate databases..."
abricate --setupdb 2>/dev/null || echo "ABRicate setup skipped"

#--- VFDB Core (~50 MB) ---
echo "Downloading VFDB core..."
if [ ! -f "${DB_ROOT}/vfdb/VFDB_setA_pro.fas" ]; then
    mkdir -p "${DB_ROOT}/vfdb"
    cd "${DB_ROOT}/vfdb"
    wget -q "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz" && gunzip -f VFDB_setA_pro.fas.gz
    wget -q "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz" && gunzip -f VFDB_setA_nt.fas.gz
    makeblastdb -in VFDB_setA_pro.fas -dbtype prot -out vfdb_prot 2>/dev/null
    makeblastdb -in VFDB_setA_nt.fas -dbtype nucl -out vfdb_nucl 2>/dev/null
else
    echo "  Already exists, skipping."
fi

#-------------------------------------------------------------------------------
# STEP 3: Create config file
#-------------------------------------------------------------------------------
echo ""
echo "=== STEP 3: Creating Configuration ==="

cat > "${DB_ROOT}/config.yaml" << EOF
# Ultra-minimal laptop configuration
# Heavy analysis (assembly, Kraken2, GTDB-Tk) should run on SERVER

databases:
  amrfinder: "${DB_ROOT}/amrfinder/latest"
  vfdb: "${DB_ROOT}/vfdb"
  # mlst uses its built-in path
  # abricate uses its built-in databases

# These analyses run LOCALLY on laptop:
local_analyses:
  - amr_detection      # AMRFinderPlus, ABRicate
  - mlst_typing        # MLST
  - virulence_factors  # VFDB
  - qc_reports         # FastQC, MultiQC, QUAST
  - report_generation  # HTML reports

# These analyses should run on SERVER (results imported):
server_analyses:
  - assembly           # SPAdes, Flye
  - taxonomy           # Kraken2 (standard), GTDB-Tk
  - full_annotation    # Prokka with full databases
EOF

echo "Config saved to: ${DB_ROOT}/config.yaml"

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
echo ""
echo "==============================================================================="
echo " SETUP COMPLETE"
echo "==============================================================================="
echo ""
echo "Disk usage:"
du -sh "${DB_ROOT}" 2>/dev/null || echo "  Could not calculate"
conda info --envs | grep gp-minimal || echo ""
echo ""
echo "To use:"
echo "  conda activate gp-minimal"
echo ""
echo "What you can do on this laptop:"
echo "  ✓ AMR gene detection (amrfinder, abricate)"
echo "  ✓ MLST typing"
echo "  ✓ Virulence factor screening"
echo "  ✓ QC reports (FastQC, MultiQC)"
echo "  ✓ Assembly QC (QUAST) - bring assembly from server"
echo "  ✓ Report generation"
echo "  ✓ Pipeline development"
echo ""
echo "What to run on SERVER:"
echo "  → Assembly (SPAdes)"
echo "  → Kraken2 classification"
echo "  → GTDB-Tk taxonomy"
echo "  → Prokka annotation"
echo ""
echo "WORKFLOW: Run assembly/taxonomy on server → Transfer results → "
echo "          Run AMR/MLST/reports on laptop"
echo ""
