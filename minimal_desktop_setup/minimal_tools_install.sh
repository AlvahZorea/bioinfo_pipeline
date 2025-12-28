#!/bin/bash
#===============================================================================
# GENOME ANALYSIS PIPELINE - MINIMAL TOOLS INSTALLATION
#===============================================================================
# Installs the essential tools needed for desktop testing via conda/mamba
# 
# Prerequisites: conda or mamba installed
#
# Usage:
#   chmod +x minimal_tools_install.sh
#   ./minimal_tools_install.sh
#===============================================================================

set -e

# Check for conda/mamba
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
elif command -v conda &> /dev/null; then
    CONDA_CMD="conda"
else
    echo "ERROR: Neither conda nor mamba found. Please install miniconda/mambaforge first."
    echo "       https://github.com/conda-forge/miniforge#mambaforge"
    exit 1
fi

echo "Using: $CONDA_CMD"
echo ""

#-------------------------------------------------------------------------------
# Option 1: Single environment (simpler but larger)
#-------------------------------------------------------------------------------
create_unified_env() {
    echo "Creating unified genome-pipeline environment..."
    
    $CONDA_CMD create -n genome-pipeline -y \
        -c conda-forge -c bioconda \
        python=3.11 \
        snakemake=8 \
        biopython \
        pandas \
        numpy \
        plotly \
        jinja2 \
        pyyaml \
        fastqc \
        multiqc \
        fastp \
        spades \
        flye \
        quast \
        samtools \
        minimap2 \
        bwa \
        kraken2 \
        bakta \
        prokka \
        ncbi-amrfinderplus \
        mlst \
        abricate \
        blast
    
    echo ""
    echo "Unified environment created: genome-pipeline"
    echo "Activate with: conda activate genome-pipeline"
}

#-------------------------------------------------------------------------------
# Option 2: Modular environments (recommended for flexibility)
#-------------------------------------------------------------------------------
create_modular_envs() {
    echo "Creating modular environments..."
    
    # Core environment (always needed)
    echo "1/5 - Creating core environment..."
    $CONDA_CMD create -n gp-core -y \
        -c conda-forge -c bioconda \
        python=3.11 \
        snakemake=8 \
        biopython \
        pandas \
        numpy \
        scipy \
        plotly \
        jinja2 \
        pyyaml \
        pytest
    
    # QC environment
    echo "2/5 - Creating QC environment..."
    $CONDA_CMD create -n gp-qc -y \
        -c conda-forge -c bioconda \
        fastqc \
        multiqc \
        fastp \
        seqkit
    
    # Assembly environment
    echo "3/5 - Creating assembly environment..."
    $CONDA_CMD create -n gp-assembly -y \
        -c conda-forge -c bioconda \
        spades \
        flye \
        quast \
        samtools \
        minimap2 \
        bwa
    
    # Annotation environment
    echo "4/5 - Creating annotation environment..."
    $CONDA_CMD create -n gp-annotation -y \
        -c conda-forge -c bioconda \
        bakta \
        prokka \
        kraken2
    
    # AMR/Typing environment
    echo "5/5 - Creating AMR/typing environment..."
    $CONDA_CMD create -n gp-amr -y \
        -c conda-forge -c bioconda \
        ncbi-amrfinderplus \
        mlst \
        abricate \
        blast
    
    echo ""
    echo "Modular environments created:"
    echo "  - gp-core      : Snakemake and Python packages"
    echo "  - gp-qc        : FastQC, MultiQC, fastp"
    echo "  - gp-assembly  : SPAdes, Flye, QUAST"
    echo "  - gp-annotation: Bakta, Prokka, Kraken2"
    echo "  - gp-amr       : AMRFinder, MLST, ABRicate"
}

#-------------------------------------------------------------------------------
# Main
#-------------------------------------------------------------------------------
echo "==============================================================================="
echo " MINIMAL TOOLS INSTALLATION"
echo "==============================================================================="
echo ""
echo "Choose installation method:"
echo "  1) Unified environment (simpler, ~5-8 GB)"
echo "  2) Modular environments (flexible, ~6-10 GB)"
echo "  3) Skip (I'll install tools manually)"
echo ""
read -p "Selection [1/2/3]: " choice

case $choice in
    1) create_unified_env ;;
    2) create_modular_envs ;;
    3) echo "Skipping tool installation." ;;
    *) echo "Invalid choice. Exiting."; exit 1 ;;
esac

echo ""
echo "==============================================================================="
echo " NEXT STEPS"
echo "==============================================================================="
echo ""
echo "1. Activate your environment:"
echo "   conda activate genome-pipeline  # (or gp-core for modular)"
echo ""
echo "2. Download minimal databases:"
echo "   ./download_minimal_databases.sh ~/genome_pipeline_db"
echo ""
echo "3. Verify installation:"
echo "   ./verify_minimal_setup.sh ~/genome_pipeline_db"
echo ""
