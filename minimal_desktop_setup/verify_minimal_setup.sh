#!/bin/bash
#===============================================================================
# GENOME ANALYSIS PIPELINE - MINIMAL SETUP VERIFICATION
#===============================================================================
# Quick verification that essential tools and databases are ready for testing
#
# Usage: ./verify_minimal_setup.sh /path/to/databases
#===============================================================================

DB_ROOT="${1:-$HOME/genome_pipeline_db}"

GREEN='\033[0;32m'
RED='\033[0;31m'
YELLOW='\033[1;33m'
NC='\033[0m'

PASS=0
FAIL=0
WARN=0

check() {
    local name="$1"
    local cmd="$2"
    if eval "$cmd" &>/dev/null; then
        echo -e "${GREEN}✓${NC} $name"
        ((PASS++))
    else
        echo -e "${RED}✗${NC} $name"
        ((FAIL++))
    fi
}

check_warn() {
    local name="$1"
    local cmd="$2"
    if eval "$cmd" &>/dev/null; then
        echo -e "${GREEN}✓${NC} $name"
        ((PASS++))
    else
        echo -e "${YELLOW}○${NC} $name (optional)"
        ((WARN++))
    fi
}

echo "==============================================================================="
echo " MINIMAL SETUP VERIFICATION"
echo "==============================================================================="
echo ""
echo "Database path: ${DB_ROOT}"
echo ""

#--- Essential Tools ---
echo "ESSENTIAL TOOLS:"
echo "----------------"
check "Python 3"         "python3 --version"
check "Snakemake"        "snakemake --version"
check "SPAdes"           "spades.py --version"
check "Flye"             "flye --version"
check "QUAST"            "quast.py --version"
check "samtools"         "samtools --version"
check "minimap2"         "minimap2 --version"
check "FastQC"           "fastqc --version"
check "fastp"            "fastp --version"
check "Kraken2"          "kraken2 --version"
check "Bakta"            "bakta --version"
check "AMRFinderPlus"    "amrfinder --version"
check "MLST"             "mlst --version"
check "BLAST+"           "blastn -version"
echo ""

#--- Optional Tools ---
echo "OPTIONAL TOOLS:"
echo "---------------"
check_warn "Prokka"      "prokka --version"
check_warn "ABRicate"    "abricate --version"
check_warn "Unicycler"   "unicycler --version"
check_warn "CheckV"      "checkv -h"
check_warn "RGI (CARD)"  "rgi main --version"
check_warn "MultiQC"     "multiqc --version"
echo ""

#--- Databases ---
echo "DATABASES:"
echo "----------"
check "Bakta DB"         "[ -d '${DB_ROOT}/bakta/db' ]"
check "Kraken2 DB"       "[ -f '${DB_ROOT}/kraken2_mini/hash.k2d' ] || [ -f '${DB_ROOT}/kraken2/standard/hash.k2d' ]"
check "AMRFinder DB"     "[ -d '${DB_ROOT}/amrfinder' ]"
check "VFDB"             "[ -f '${DB_ROOT}/vfdb/VFDB_setA_pro.fas' ]"
check_warn "CheckV DB"   "[ -d '${DB_ROOT}/checkv' ]"
check_warn "CARD DB"     "[ -f '${DB_ROOT}/card/card.json' ]"
echo ""

#--- Python Packages ---
echo "PYTHON PACKAGES:"
echo "----------------"
check "BioPython"        "python3 -c 'import Bio'"
check "pandas"           "python3 -c 'import pandas'"
check "numpy"            "python3 -c 'import numpy'"
check "plotly"           "python3 -c 'import plotly'"
check "jinja2"           "python3 -c 'import jinja2'"
check "pyyaml"           "python3 -c 'import yaml'"
echo ""

#--- Summary ---
echo "==============================================================================="
echo " SUMMARY"
echo "==============================================================================="
echo ""
echo -e "Passed:   ${GREEN}${PASS}${NC}"
echo -e "Failed:   ${RED}${FAIL}${NC}"
echo -e "Optional: ${YELLOW}${WARN}${NC} (not required for basic testing)"
echo ""

if [ $FAIL -eq 0 ]; then
    echo -e "${GREEN}Ready for testing!${NC}"
    echo ""
    echo "You can now:"
    echo "  1. Run the pipeline on test data"
    echo "  2. Test individual modules"
    echo ""
else
    echo -e "${RED}Some required components are missing.${NC}"
    echo ""
    echo "To fix:"
    echo "  1. Install missing tools via conda"
    echo "  2. Download missing databases"
    echo "  3. Re-run this verification"
fi
echo ""
