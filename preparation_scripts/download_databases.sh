#!/bin/bash
#===============================================================================
# GENOME ANALYSIS PIPELINE - DATABASE DOWNLOAD SCRIPT
#===============================================================================
# Run this script on your server this week to download all required databases
# Estimated total download size: ~200-350 GB
# Estimated time: Several hours depending on connection speed
#
# Usage: 
#   chmod +x download_databases.sh
#   ./download_databases.sh /path/to/database/directory
#===============================================================================

set -e  # Exit on error

# Configuration
DB_ROOT="${1:-/databases}"
THREADS="${2:-8}"

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

echo_status() {
    echo -e "${GREEN}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} $1"
}

echo_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

echo_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

#-------------------------------------------------------------------------------
# Create directory structure
#-------------------------------------------------------------------------------
echo_status "Creating database directory structure..."

mkdir -p "${DB_ROOT}"/{kraken2,gtdbtk,bakta,checkv,genomad,card,amrfinder}
mkdir -p "${DB_ROOT}"/{vfdb,pubmlst,pharokka,interpro,eggnog,phrogs,iceberg}
mkdir -p "${DB_ROOT}"/{centrifuge,refseq_viruses,ncbi_taxonomy}

#-------------------------------------------------------------------------------
# PRIORITY 1: Essential Databases
#-------------------------------------------------------------------------------

echo ""
echo "==============================================================================="
echo " PRIORITY 1: ESSENTIAL DATABASES (~200 GB)"
echo "==============================================================================="
echo ""

#--- 1. Kraken2 Standard Database (~50-100 GB) ---
echo_status "1/10 - Downloading Kraken2 Standard Database..."
echo "       This is the largest download and may take several hours."
echo "       Size: ~50-100 GB"

if [ ! -f "${DB_ROOT}/kraken2/standard/hash.k2d" ]; then
    cd "${DB_ROOT}/kraken2"
    # Option A: Download pre-built (faster)
    # wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20230605.tar.gz
    # tar -xzf k2_standard_20230605.tar.gz -C standard/
    
    # Option B: Build from scratch (more current)
    kraken2-build --standard --threads ${THREADS} --db "${DB_ROOT}/kraken2/standard" || \
        echo_warning "Kraken2 build failed - may need more memory or disk space"
    echo_status "Kraken2 database complete."
else
    echo_status "Kraken2 database already exists, skipping."
fi

#--- 2. Bakta Database (~30 GB) ---
echo_status "2/10 - Downloading Bakta Database..."
echo "       Size: ~30 GB (full) or ~2 GB (light)"

if [ ! -d "${DB_ROOT}/bakta/db" ]; then
    # Full version (recommended)
    bakta_db download --output "${DB_ROOT}/bakta" --type full || \
        echo_error "Bakta download failed - check bakta installation"
    echo_status "Bakta database complete."
else
    echo_status "Bakta database already exists, skipping."
fi

#--- 3. GTDB-Tk Database (~85 GB) ---
echo_status "3/10 - Downloading GTDB-Tk Database..."
echo "       Size: ~85 GB"

if [ ! -d "${DB_ROOT}/gtdbtk/release" ]; then
    cd "${DB_ROOT}/gtdbtk"
    # Download latest release
    wget -c https://data.gtdb.ecogenomic.org/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz || \
        echo_error "GTDB-Tk download failed - check URL for latest version"
    tar -xzf gtdbtk_data.tar.gz
    rm gtdbtk_data.tar.gz
    echo_status "GTDB-Tk database complete."
else
    echo_status "GTDB-Tk database already exists, skipping."
fi

#--- 4. CheckV Database (~2 GB) ---
echo_status "4/10 - Downloading CheckV Database..."
echo "       Size: ~2 GB"

if [ ! -d "${DB_ROOT}/checkv/checkv-db" ]; then
    checkv download_database "${DB_ROOT}/checkv" || \
        echo_error "CheckV download failed - check checkv installation"
    echo_status "CheckV database complete."
else
    echo_status "CheckV database already exists, skipping."
fi

#--- 5. geNomad Database (~3 GB) ---
echo_status "5/10 - Downloading geNomad Database..."
echo "       Size: ~3 GB"

if [ ! -d "${DB_ROOT}/genomad/genomad_db" ]; then
    genomad download-database "${DB_ROOT}/genomad" || \
        echo_error "geNomad download failed - check genomad installation"
    echo_status "geNomad database complete."
else
    echo_status "geNomad database already exists, skipping."
fi

#--- 6. CARD Database (~500 MB) ---
echo_status "6/10 - Downloading CARD Database..."
echo "       Size: ~500 MB"

if [ ! -f "${DB_ROOT}/card/card.json" ]; then
    cd "${DB_ROOT}/card"
    wget -c https://card.mcmaster.ca/latest/data -O card-data.tar.bz2
    tar -xjf card-data.tar.bz2
    rm card-data.tar.bz2
    echo_status "CARD database complete."
else
    echo_status "CARD database already exists, skipping."
fi

#--- 7. AMRFinderPlus Database (~500 MB) ---
echo_status "7/10 - Downloading AMRFinderPlus Database..."
echo "       Size: ~500 MB"

if [ ! -d "${DB_ROOT}/amrfinder/latest" ]; then
    amrfinder_update -d "${DB_ROOT}/amrfinder" || \
        echo_error "AMRFinderPlus update failed - check installation"
    echo_status "AMRFinderPlus database complete."
else
    echo_status "AMRFinderPlus database already exists, skipping."
fi

#--- 8. VFDB Database (~100 MB) ---
echo_status "8/10 - Downloading VFDB Database..."
echo "       Size: ~100 MB"

if [ ! -f "${DB_ROOT}/vfdb/VFDB_setA_pro.fas" ]; then
    cd "${DB_ROOT}/vfdb"
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setA_nt.fas.gz
    wget -c http://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
    gunzip *.gz
    
    # Create BLAST databases
    makeblastdb -in VFDB_setA_pro.fas -dbtype prot -out vfdb_core_prot -title "VFDB Core Proteins"
    makeblastdb -in VFDB_setB_pro.fas -dbtype prot -out vfdb_full_prot -title "VFDB Full Proteins"
    makeblastdb -in VFDB_setA_nt.fas -dbtype nucl -out vfdb_core_nucl -title "VFDB Core Nucleotides"
    makeblastdb -in VFDB_setB_nt.fas -dbtype nucl -out vfdb_full_nucl -title "VFDB Full Nucleotides"
    echo_status "VFDB database complete."
else
    echo_status "VFDB database already exists, skipping."
fi

#--- 9. PubMLST Schemes (~1 GB) ---
echo_status "9/10 - Downloading PubMLST MLST Schemes..."
echo "       Size: ~1 GB"

# mlst tool handles its own database
mlst --update 2>/dev/null || echo_warning "mlst update failed - may need manual setup"
echo_status "PubMLST schemes update complete."

#--- 10. Pharokka Database (~8 GB) ---
echo_status "10/10 - Downloading Pharokka Database..."
echo "        Size: ~8 GB"

if [ ! -d "${DB_ROOT}/pharokka/pharokka_db" ]; then
    # Pharokka's install script
    install_databases.py -o "${DB_ROOT}/pharokka" 2>/dev/null || \
        echo_warning "Pharokka database download failed - check pharokka installation"
    echo_status "Pharokka database complete."
else
    echo_status "Pharokka database already exists, skipping."
fi

#-------------------------------------------------------------------------------
# PRIORITY 2: Recommended Databases (Optional but useful)
#-------------------------------------------------------------------------------

echo ""
echo "==============================================================================="
echo " PRIORITY 2: RECOMMENDED DATABASES (~150 GB) - Optional"
echo "==============================================================================="
echo ""
read -p "Download Priority 2 databases? (y/n): " download_p2

if [ "$download_p2" == "y" ]; then

    #--- InterProScan Data (~50 GB) ---
    echo_status "Downloading InterProScan Data..."
    echo "       Size: ~50 GB (very large)"
    echo "       This is optional but highly useful for domain annotation"
    
    if [ ! -d "${DB_ROOT}/interpro/data" ]; then
        cd "${DB_ROOT}/interpro"
        # Check InterPro website for latest version URL
        INTERPRO_VERSION="5.66-98.0"
        wget -c "https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${INTERPRO_VERSION}/interproscan-${INTERPRO_VERSION}-64-bit.tar.gz"
        tar -xzf "interproscan-${INTERPRO_VERSION}-64-bit.tar.gz"
        echo_status "InterProScan complete."
    else
        echo_status "InterProScan already exists, skipping."
    fi

    #--- eggNOG Database (~40 GB) ---
    echo_status "Downloading eggNOG Database..."
    echo "       Size: ~40 GB"
    
    if [ ! -d "${DB_ROOT}/eggnog/eggnog.db" ]; then
        download_eggnog_data.py -y --data_dir "${DB_ROOT}/eggnog" || \
            echo_warning "eggNOG download failed"
        echo_status "eggNOG database complete."
    else
        echo_status "eggNOG database already exists, skipping."
    fi

    #--- PHROGs Database (~1 GB) ---
    echo_status "Downloading PHROGs Database..."
    echo "       Size: ~1 GB"
    
    if [ ! -d "${DB_ROOT}/phrogs" ]; then
        cd "${DB_ROOT}/phrogs"
        wget -c https://phrogs.lmge.uca.fr/downloads/phrogs_hhsuite_db.tar.gz
        tar -xzf phrogs_hhsuite_db.tar.gz
        echo_status "PHROGs database complete."
    else
        echo_status "PHROGs database already exists, skipping."
    fi

fi

#-------------------------------------------------------------------------------
# Create configuration file with database paths
#-------------------------------------------------------------------------------

echo_status "Creating database paths configuration file..."

cat > "${DB_ROOT}/database_paths.yaml" << EOF
# Database paths configuration
# Generated on: $(date)
# Edit paths as needed for your installation

databases:
  kraken2_standard: "${DB_ROOT}/kraken2/standard"
  gtdbtk: "${DB_ROOT}/gtdbtk/release*/release*"
  bakta: "${DB_ROOT}/bakta/db"
  checkv: "${DB_ROOT}/checkv/checkv-db-v*"
  genomad: "${DB_ROOT}/genomad/genomad_db"
  card: "${DB_ROOT}/card"
  amrfinder: "${DB_ROOT}/amrfinder/latest"
  vfdb: "${DB_ROOT}/vfdb"
  pharokka: "${DB_ROOT}/pharokka"
  interpro: "${DB_ROOT}/interpro/interproscan-*"
  eggnog: "${DB_ROOT}/eggnog"
  phrogs: "${DB_ROOT}/phrogs"
EOF

echo_status "Configuration file created: ${DB_ROOT}/database_paths.yaml"

#-------------------------------------------------------------------------------
# Validation
#-------------------------------------------------------------------------------

echo ""
echo "==============================================================================="
echo " DATABASE INSTALLATION SUMMARY"
echo "==============================================================================="
echo ""

validate_db() {
    local name="$1"
    local path="$2"
    if [ -e "$path" ]; then
        echo -e "${GREEN}✓${NC} $name: OK"
        return 0
    else
        echo -e "${RED}✗${NC} $name: MISSING"
        return 1
    fi
}

validate_db "Kraken2 Standard" "${DB_ROOT}/kraken2/standard/hash.k2d"
validate_db "Bakta" "${DB_ROOT}/bakta/db"
validate_db "GTDB-Tk" "${DB_ROOT}/gtdbtk/release"*
validate_db "CheckV" "${DB_ROOT}/checkv/checkv-db"*
validate_db "geNomad" "${DB_ROOT}/genomad/genomad_db"
validate_db "CARD" "${DB_ROOT}/card/card.json"
validate_db "AMRFinderPlus" "${DB_ROOT}/amrfinder/latest"
validate_db "VFDB" "${DB_ROOT}/vfdb/VFDB_setA_pro.fas"
validate_db "Pharokka" "${DB_ROOT}/pharokka"

echo ""
echo_status "Database download script complete!"
echo ""
echo "Total disk usage:"
du -sh "${DB_ROOT}"/* 2>/dev/null || echo "Could not calculate disk usage"
echo ""
echo "Next steps:"
echo "1. Verify all databases are correctly downloaded"
echo "2. Note any failures above and retry manually if needed"
echo "3. Run validation tests (see validate_installation.sh)"
echo "4. Transfer databases to air-gapped server if needed"
echo ""

