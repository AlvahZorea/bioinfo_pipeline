#!/bin/bash
#===============================================================================
# GENOME ANALYSIS PIPELINE - MINIMAL DESKTOP TESTING DATABASE SETUP
#===============================================================================
# This is a LIGHTWEIGHT version for local testing on a desktop/laptop
# Total size: ~15-25 GB (vs ~200-350 GB for full production)
#
# Usage: 
#   chmod +x download_minimal_databases.sh
#   ./download_minimal_databases.sh /path/to/databases
#
# Note: This setup is for TESTING ONLY. For production/air-gapped deployment,
#       use the full download_databases.sh script.
#===============================================================================

set -e

# Configuration
DB_ROOT="${1:-$HOME/genome_pipeline_db}"
THREADS="${2:-4}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

echo_status() { echo -e "${GREEN}[INFO]${NC} $1"; }
echo_warning() { echo -e "${YELLOW}[WARN]${NC} $1"; }
echo_error() { echo -e "${RED}[ERROR]${NC} $1"; }
echo_section() { echo -e "\n${CYAN}=== $1 ===${NC}\n"; }

#-------------------------------------------------------------------------------
# Create directory structure
#-------------------------------------------------------------------------------
echo_section "Setting up minimal database directory"

mkdir -p "${DB_ROOT}"/{bakta,kraken2_mini,amrfinder,vfdb,mlst}
mkdir -p "${DB_ROOT}"/{checkv,card}

echo "Database root: ${DB_ROOT}"
echo ""

#-------------------------------------------------------------------------------
# TIER 1: Core Essentials (~10-15 GB)
# These are needed for basic bacterial genome analysis
#-------------------------------------------------------------------------------
echo_section "TIER 1: Core Essentials (~10-15 GB)"

#--- Bakta Light Database (~2 GB) ---
# Modern bacterial annotation - HIGHLY RECOMMENDED
echo_status "1/5 - Bakta Light Database (~2 GB)"
echo "       Used for: Gene prediction and annotation"

if [ ! -d "${DB_ROOT}/bakta/db" ]; then
    if command -v bakta_db &> /dev/null; then
        bakta_db download --output "${DB_ROOT}/bakta" --type light && \
            echo_status "Bakta Light database complete." || \
            echo_warning "Bakta download failed - will need manual setup"
    else
        echo_warning "bakta_db command not found - install bakta first or skip"
    fi
else
    echo_status "Bakta database already exists, skipping."
fi

#--- Kraken2 Mini Database (~8 GB) ---
# Taxonomic classification - use mini version for testing
echo_status "2/5 - Kraken2 Mini Database (~8 GB)"
echo "       Used for: Species identification and contamination screening"
echo "       Note: Mini DB is sufficient for common organisms in testing"

if [ ! -f "${DB_ROOT}/kraken2_mini/hash.k2d" ]; then
    cd "${DB_ROOT}/kraken2_mini"
    
    # Option 1: Download pre-built MiniKraken (smaller, older)
    # wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v2_8GB_201904.tgz
    # tar -xzf minikraken2_v2_8GB_201904.tgz
    
    # Option 2: Download k2_standard_08gb (recommended - more recent)
    echo "       Downloading pre-built 8GB database..."
    wget -c "https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20231009.tar.gz" -O k2_mini.tar.gz && \
    tar -xzf k2_mini.tar.gz && \
    rm k2_mini.tar.gz && \
    echo_status "Kraken2 mini database complete." || \
    echo_warning "Kraken2 download failed - check URL or network"
else
    echo_status "Kraken2 mini database already exists, skipping."
fi

#--- AMRFinderPlus Database (~500 MB) ---
# AMR detection - essential for clinical relevance
echo_status "3/5 - AMRFinderPlus Database (~500 MB)"
echo "       Used for: Antimicrobial resistance gene detection"

if [ ! -d "${DB_ROOT}/amrfinder/latest" ]; then
    if command -v amrfinder_update &> /dev/null; then
        amrfinder_update -d "${DB_ROOT}/amrfinder" && \
            echo_status "AMRFinderPlus database complete." || \
            echo_warning "AMRFinderPlus update failed"
    else
        echo_warning "amrfinder_update not found - install ncbi-amrfinderplus first"
    fi
else
    echo_status "AMRFinderPlus database already exists, skipping."
fi

#--- MLST Schemes (~500 MB - 1 GB) ---
# Strain typing
echo_status "4/5 - MLST Schemes (~500 MB - 1 GB)"
echo "       Used for: Multi-locus sequence typing"

if command -v mlst &> /dev/null; then
    mlst --check 2>/dev/null || mlst --update 2>/dev/null && \
        echo_status "MLST schemes updated." || \
        echo_warning "MLST update failed - may need manual setup"
else
    echo_warning "mlst command not found - install mlst first"
fi

#--- VFDB Core (~50 MB) ---
# Virulence factors - core set only
echo_status "5/5 - VFDB Core Database (~50 MB)"
echo "       Used for: Virulence factor detection"

if [ ! -f "${DB_ROOT}/vfdb/VFDB_setA_pro.fas" ]; then
    cd "${DB_ROOT}/vfdb"
    wget -c "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz" && \
    wget -c "http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz" && \
    gunzip -f *.gz && \
    echo_status "VFDB core database downloaded." || \
    echo_warning "VFDB download failed"
    
    # Create BLAST databases if makeblastdb available
    if command -v makeblastdb &> /dev/null; then
        makeblastdb -in VFDB_setA_pro.fas -dbtype prot -out vfdb_prot 2>/dev/null
        makeblastdb -in VFDB_setA_nt.fas -dbtype nucl -out vfdb_nucl 2>/dev/null
        echo_status "VFDB BLAST databases created."
    fi
else
    echo_status "VFDB database already exists, skipping."
fi

#-------------------------------------------------------------------------------
# TIER 2: Optional Add-ons (~5-10 GB additional)
# Uncomment if you have space and want more complete testing
#-------------------------------------------------------------------------------
echo_section "TIER 2: Optional Databases"

read -p "Download optional databases (CheckV for viruses, CARD for AMR)? [y/N]: " download_optional

if [[ "$download_optional" =~ ^[Yy]$ ]]; then
    
    #--- CheckV Database (~2 GB) ---
    echo_status "Optional: CheckV Database (~2 GB)"
    echo "          Used for: Viral genome quality assessment"
    
    if [ ! -d "${DB_ROOT}/checkv/checkv-db"* ]; then
        if command -v checkv &> /dev/null; then
            checkv download_database "${DB_ROOT}/checkv" && \
                echo_status "CheckV database complete." || \
                echo_warning "CheckV download failed"
        else
            echo_warning "checkv not found - skip or install first"
        fi
    fi
    
    #--- CARD Database (~500 MB) ---
    echo_status "Optional: CARD Database (~500 MB)"
    echo "          Used for: Comprehensive AMR analysis with RGI"
    
    if [ ! -f "${DB_ROOT}/card/card.json" ]; then
        cd "${DB_ROOT}/card"
        wget -c "https://card.mcmaster.ca/latest/data" -O card-data.tar.bz2 && \
        tar -xjf card-data.tar.bz2 && \
        rm card-data.tar.bz2 && \
        echo_status "CARD database complete." || \
        echo_warning "CARD download failed"
    fi
fi

#-------------------------------------------------------------------------------
# Create configuration file
#-------------------------------------------------------------------------------
echo_section "Creating Configuration"

cat > "${DB_ROOT}/database_paths.yaml" << EOF
# Minimal Desktop Testing Database Configuration
# Generated: $(date)
# 
# This is a LIGHTWEIGHT configuration for testing purposes.
# For production use, see the full database setup.

databases:
  # Tier 1: Core (required for basic testing)
  bakta: "${DB_ROOT}/bakta/db"
  kraken2: "${DB_ROOT}/kraken2_mini"
  amrfinder: "${DB_ROOT}/amrfinder/latest"
  mlst: "default"  # Uses mlst's built-in path
  vfdb: "${DB_ROOT}/vfdb"
  
  # Tier 2: Optional
  checkv: "${DB_ROOT}/checkv"
  card: "${DB_ROOT}/card"

# Note: The following are NOT included in minimal setup
# For full analysis, you'll need:
#   - gtdbtk (~85 GB) - detailed bacterial taxonomy
#   - kraken2_standard (~100 GB) - comprehensive taxonomy
#   - pharokka (~8 GB) - phage annotation
#   - genomad (~3 GB) - virus/plasmid detection
#   - interproscan (~50 GB) - protein domains
EOF

echo_status "Configuration saved to: ${DB_ROOT}/database_paths.yaml"

#-------------------------------------------------------------------------------
# Summary
#-------------------------------------------------------------------------------
echo_section "INSTALLATION SUMMARY"

echo "Database location: ${DB_ROOT}"
echo ""
echo "Disk usage:"
du -sh "${DB_ROOT}"/* 2>/dev/null | sort -h || echo "  (could not calculate)"
echo ""
echo "Total:"
du -sh "${DB_ROOT}" 2>/dev/null || echo "  (could not calculate)"
echo ""

echo_status "Minimal desktop setup complete!"
echo ""
echo "What you can test with this setup:"
echo "  ✓ Genome assembly (SPAdes, Flye, Unicycler)"
echo "  ✓ Basic taxonomic ID (Kraken2 mini)"
echo "  ✓ Gene annotation (Bakta light)"
echo "  ✓ AMR detection (AMRFinderPlus)"
echo "  ✓ MLST typing"
echo "  ✓ Virulence factors (VFDB core)"
echo ""
echo "What requires full databases (run on server):"
echo "  ⚠ Detailed taxonomy (GTDB-Tk)"
echo "  ⚠ Rare organism identification"
echo "  ⚠ Phage annotation (Pharokka)"
echo "  ⚠ Comprehensive protein domains (InterProScan)"
echo ""
