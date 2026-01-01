#!/bin/bash
#===============================================================================
# GENOME ANALYSIS PIPELINE - INSTALLATION VALIDATION SCRIPT
#===============================================================================
# Run this script to verify all tools and databases are correctly installed
# before transferring to the air-gapped server
#
# Usage: ./validate_installation.sh /path/to/database/directory
#===============================================================================

DB_ROOT="${1:-/databases}"
TEST_DIR="${2:-/tmp/pipeline_validation}"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

PASS_COUNT=0
FAIL_COUNT=0
WARN_COUNT=0

check_pass() {
    echo -e "${GREEN}✓ PASS${NC}: $1"
    ((PASS_COUNT++))
}

check_fail() {
    echo -e "${RED}✗ FAIL${NC}: $1"
    ((FAIL_COUNT++))
}

check_warn() {
    echo -e "${YELLOW}⚠ WARN${NC}: $1"
    ((WARN_COUNT++))
}

check_cmd() {
    local cmd="$1"
    local name="${2:-$cmd}"
    if command -v "$cmd" &> /dev/null; then
        local version=$($cmd --version 2>&1 | head -1 || echo "version unknown")
        check_pass "$name ($version)"
    else
        check_fail "$name not found"
    fi
}

check_db() {
    local name="$1"
    local path="$2"
    if [ -e "$path" ]; then
        local size=$(du -sh "$path" 2>/dev/null | cut -f1)
        check_pass "$name (${size})"
    else
        check_fail "$name not found at $path"
    fi
}

echo "==============================================================================="
echo " GENOME ANALYSIS PIPELINE - INSTALLATION VALIDATION"
echo "==============================================================================="
echo ""
echo "Database root: ${DB_ROOT}"
echo "Test directory: ${TEST_DIR}"
echo ""

#-------------------------------------------------------------------------------
# 1. Check Core Tools
#-------------------------------------------------------------------------------
echo -e "${BLUE}[1/6] Checking Core Tools${NC}"
echo "----------------------------------------"

# Quality Control
check_cmd "fastqc" "FastQC"
check_cmd "multiqc" "MultiQC"
check_cmd "fastp" "fastp"

# Assembly
check_cmd "spades.py" "SPAdes"
check_cmd "flye" "Flye"
check_cmd "unicycler" "Unicycler"
check_cmd "medaka" "Medaka"
check_cmd "quast.py" "QUAST"

# Mapping
check_cmd "minimap2" "minimap2"
check_cmd "bwa" "BWA"
check_cmd "samtools" "samtools"

# Identification
check_cmd "kraken2" "Kraken2"
check_cmd "gtdbtk" "GTDB-Tk"
check_cmd "mash" "Mash"
check_cmd "fastANI" "fastANI"
check_cmd "genomad" "geNomad"

# Annotation
check_cmd "prokka" "Prokka"
check_cmd "bakta" "Bakta"
check_cmd "pharokka" "Pharokka"

# AMR
check_cmd "amrfinder" "AMRFinderPlus"
check_cmd "abricate" "ABRicate"
check_cmd "rgi" "RGI (CARD)"

# Typing
check_cmd "mlst" "MLST"

# MGE Detection
check_cmd "mob_recon" "MOB-suite"
check_cmd "phispy" "PhiSpy"

# Phylogeny
check_cmd "mafft" "MAFFT"
check_cmd "iqtree2" "IQ-TREE2"
check_cmd "parsnp" "Parsnp"

# Utilities
check_cmd "blastn" "BLAST+"
check_cmd "makeblastdb" "BLAST+ (makeblastdb)"
check_cmd "snakemake" "Snakemake"
check_cmd "python3" "Python3"

echo ""

#-------------------------------------------------------------------------------
# 2. Check Databases
#-------------------------------------------------------------------------------
echo -e "${BLUE}[2/6] Checking Databases${NC}"
echo "----------------------------------------"

check_db "Kraken2 Standard" "${DB_ROOT}/kraken2/standard/hash.k2d"
check_db "Bakta" "${DB_ROOT}/bakta/db"
check_db "GTDB-Tk" "${DB_ROOT}/gtdbtk"
check_db "CheckV" "${DB_ROOT}/checkv"
check_db "geNomad" "${DB_ROOT}/genomad"
check_db "CARD" "${DB_ROOT}/card/card.json"
check_db "AMRFinderPlus" "${DB_ROOT}/amrfinder"
check_db "VFDB" "${DB_ROOT}/vfdb/VFDB_setA_pro.fas"
check_db "Pharokka" "${DB_ROOT}/pharokka"

# Optional databases
echo ""
echo "Optional databases:"
[ -e "${DB_ROOT}/interpro" ] && check_pass "InterProScan" || check_warn "InterProScan (optional)"
[ -e "${DB_ROOT}/eggnog" ] && check_pass "eggNOG" || check_warn "eggNOG (optional)"

echo ""

#-------------------------------------------------------------------------------
# 3. Run Minimal Functional Tests
#-------------------------------------------------------------------------------
echo -e "${BLUE}[3/6] Running Functional Tests${NC}"
echo "----------------------------------------"

mkdir -p "${TEST_DIR}"
cd "${TEST_DIR}"

# Create tiny test sequence
echo ">test_seq
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
GCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCAT" > test.fasta

# Create tiny FASTQ
echo "@read1
ATGCATGCATGCATGCATGC
+
FFFFFFFFFFFFFFFFFFFF
@read2
GCATGCATGCATGCATGCAT
+
FFFFFFFFFFFFFFFFFFFF" > test.fastq

# Test blastn
if blastn -version &>/dev/null; then
    makeblastdb -in test.fasta -dbtype nucl -out testdb &>/dev/null && \
    blastn -query test.fasta -db testdb -outfmt 6 &>/dev/null && \
    check_pass "BLAST+ functional test" || check_fail "BLAST+ functional test"
fi

# Test Kraken2 (requires database)
if [ -e "${DB_ROOT}/kraken2/standard/hash.k2d" ]; then
    kraken2 --db "${DB_ROOT}/kraken2/standard" --threads 1 test.fastq &>/dev/null && \
    check_pass "Kraken2 functional test" || check_fail "Kraken2 functional test"
else
    check_warn "Skipping Kraken2 test (no database)"
fi

# Test Python imports
python3 -c "import Bio; import pandas; import numpy; import plotly; import snakemake" 2>/dev/null && \
    check_pass "Python core packages" || check_fail "Python core packages"

# Cleanup
rm -rf "${TEST_DIR}"

echo ""

#-------------------------------------------------------------------------------
# 4. Check Docker Images (if using Docker)
#-------------------------------------------------------------------------------
echo -e "${BLUE}[4/6] Checking Docker Images${NC}"
echo "----------------------------------------"

if command -v docker &> /dev/null; then
    check_pass "Docker installed"
    
    # Check for key images
    docker image inspect staphb/spades:3.15.5 &>/dev/null && \
        check_pass "SPAdes Docker image" || check_warn "SPAdes Docker image not found"
    docker image inspect staphb/kraken2:2.1.3 &>/dev/null && \
        check_pass "Kraken2 Docker image" || check_warn "Kraken2 Docker image not found"
else
    check_warn "Docker not installed (optional)"
fi

echo ""

#-------------------------------------------------------------------------------
# 5. Check Disk Space
#-------------------------------------------------------------------------------
echo -e "${BLUE}[5/6] Checking Disk Space${NC}"
echo "----------------------------------------"

DB_DISK_FREE=$(df -BG "${DB_ROOT}" 2>/dev/null | tail -1 | awk '{print $4}' | tr -d 'G')
if [ -n "$DB_DISK_FREE" ] && [ "$DB_DISK_FREE" -gt 50 ]; then
    check_pass "Database disk space (${DB_DISK_FREE}G free)"
else
    check_warn "Low disk space on database partition (${DB_DISK_FREE}G free)"
fi

echo ""

#-------------------------------------------------------------------------------
# 6. Summary
#-------------------------------------------------------------------------------
echo -e "${BLUE}[6/6] Summary${NC}"
echo "----------------------------------------"
echo ""
echo -e "Tests passed: ${GREEN}${PASS_COUNT}${NC}"
echo -e "Tests failed: ${RED}${FAIL_COUNT}${NC}"
echo -e "Warnings:     ${YELLOW}${WARN_COUNT}${NC}"
echo ""

if [ $FAIL_COUNT -eq 0 ]; then
    echo -e "${GREEN}All critical checks passed!${NC}"
    echo "Your installation is ready for pipeline development."
else
    echo -e "${RED}Some critical checks failed.${NC}"
    echo "Please resolve the issues above before proceeding."
fi

echo ""
echo "==============================================================================="

