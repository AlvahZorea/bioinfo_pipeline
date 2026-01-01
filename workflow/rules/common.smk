"""
============================================================================
Common Rules and Utility Functions
============================================================================
Shared functions for the genome analysis pipeline
"""

import os
import sys
from pathlib import Path

# ============================================================================
# Helper Functions
# ============================================================================

def get_sample_id(config):
    """Get sample ID from config"""
    return config.get("sample_id", "sample")


def get_output_dir(config):
    """Get base output directory from config"""
    return config.get("output_dir", "output")


def get_sample_output(config):
    """Get sample-specific output directory"""
    return os.path.join(get_output_dir(config), get_sample_id(config))


def get_threads(config):
    """Get number of threads from config"""
    return config.get("threads", 8)


def get_memory_gb(config):
    """Get memory limit in GB from config"""
    return config.get("memory_gb", 16)


def get_short_reads(config):
    """
    Get short read paths from config.
    Returns tuple (r1_path, r2_path) or (None, None) if not provided.
    """
    reads = config.get("reads", {})
    short = reads.get("short", {})
    r1 = short.get("r1")
    r2 = short.get("r2")
    return r1, r2


def get_long_reads(config):
    """
    Get long read path from config.
    Returns path or None if not provided.
    """
    reads = config.get("reads", {})
    return reads.get("long")


def has_short_reads(config):
    """Check if short reads are provided"""
    r1, r2 = get_short_reads(config)
    return r1 is not None and r2 is not None


def has_long_reads(config):
    """Check if long reads are provided"""
    return get_long_reads(config) is not None


def detect_read_mode(config):
    """
    Auto-detect execution mode based on available reads.
    Returns: "short", "long", "hybrid", or raises error
    """
    mode = config.get("mode", "auto")
    
    if mode != "auto":
        return mode
    
    short = has_short_reads(config)
    long = has_long_reads(config)
    
    if short and long:
        return "hybrid"
    elif long:
        return "long"
    elif short:
        return "short"
    else:
        raise ValueError("No reads provided in configuration. Set reads.short.r1/r2 or reads.long")


def is_module_enabled(config, module_name):
    """Check if a module is enabled in config"""
    modules = config.get("modules", {})
    return modules.get(module_name, True)


def get_database_path(config, db_name, subkey=None):
    """
    Get database path from config.
    
    Args:
        config: Configuration dictionary
        db_name: Database name (e.g., "kraken2", "bakta")
        subkey: Optional subkey for nested databases (e.g., "standard" for kraken2)
    
    Returns:
        Path string or None if not found
    """
    # First try config databases
    db_config = config.get("databases", {})
    
    # Try to load from databases.yaml if not in main config
    if not db_config:
        db_file = "config/databases.yaml"
        if os.path.exists(db_file):
            import yaml
            with open(db_file) as f:
                db_data = yaml.safe_load(f)
                db_config = db_data.get("databases", {})
    
    # Get the database entry
    db_entry = db_config.get(db_name)
    
    if db_entry is None:
        return None
    
    # Handle nested databases (e.g., kraken2.standard)
    if subkey and isinstance(db_entry, dict):
        return db_entry.get(subkey)
    elif isinstance(db_entry, dict) and not subkey:
        # Return first value if it's a dict and no subkey specified
        return next(iter(db_entry.values()), None)
    else:
        return db_entry


def get_qc_params(config):
    """Get QC parameters with defaults"""
    qc = config.get("qc", {})
    return {
        "min_quality": qc.get("min_quality", 20),
        "min_length": qc.get("min_length", 50),
        "adapter_trimming": qc.get("adapter_trimming", True),
        "contamination_screening": qc.get("contamination_screening", True),
        "remove_contaminants": qc.get("remove_contaminants", True),
        "contaminant_db": qc.get("contaminant_db"),
        "nanopore": {
            "min_length": qc.get("nanopore", {}).get("min_length", 500),
            "min_quality": qc.get("nanopore", {}).get("min_quality", 7),
        }
    }


def get_assembly_params(config):
    """Get assembly parameters with defaults"""
    assembly = config.get("assembly", {})
    return {
        "assembler": assembly.get("assembler", "auto"),
        "min_contig_length": assembly.get("min_contig_length", 500),
        "polish": assembly.get("polish", True),
        "polish_rounds": assembly.get("polish_rounds", 1),
        "spades": {
            "mode": assembly.get("spades", {}).get("mode", "careful"),
            "kmer_sizes": assembly.get("spades", {}).get("kmer_sizes", [21, 33, 55, 77]),
        },
        "flye": {
            "read_type": assembly.get("flye", {}).get("read_type", "nano-raw"),
            "genome_size": assembly.get("flye", {}).get("genome_size", "5m"),
        }
    }


def select_assembler(config):
    """
    Select appropriate assembler based on read mode and config.
    Returns: "spades", "flye", or "unicycler"
    """
    assembly_params = get_assembly_params(config)
    assembler = assembly_params["assembler"]
    
    if assembler != "auto":
        return assembler
    
    mode = detect_read_mode(config)
    
    if mode == "short":
        return "spades"
    elif mode == "long":
        return "flye"
    elif mode == "hybrid":
        return "unicycler"
    else:
        return "spades"


# ============================================================================
# Output Path Functions
# ============================================================================

def qc_output(config, filename):
    """Get path for QC output file"""
    return os.path.join(get_sample_output(config), "01_qc", filename)


def assembly_output(config, filename):
    """Get path for assembly output file"""
    return os.path.join(get_sample_output(config), "02_assembly", filename)


def identification_output(config, filename):
    """Get path for identification output file"""
    return os.path.join(get_sample_output(config), "03_identification", filename)


def annotation_output(config, filename):
    """Get path for annotation output file"""
    return os.path.join(get_sample_output(config), "04_annotation", filename)


def specialized_output(config, subdir, filename):
    """Get path for specialized analysis output file"""
    return os.path.join(get_sample_output(config), "05_specialized", subdir, filename)


def report_output(config, filename):
    """Get path for report output file"""
    return os.path.join(get_sample_output(config), "08_report", filename)


def anomaly_output(config, filename):
    """Get path for anomaly detection output file"""
    return os.path.join(get_sample_output(config), "06_anomaly", filename)


def log_output(config, filename):
    """Get path for log file"""
    return os.path.join(get_sample_output(config), "00_logs", filename)
