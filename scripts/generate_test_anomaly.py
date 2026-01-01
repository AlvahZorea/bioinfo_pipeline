#!/usr/bin/env python3
"""
Generate test anomaly detection data for the Citrobacter sample.
This creates mock anomaly data to test the visualization without running BLAST.
"""

import json
import os
from pathlib import Path

# Configuration
SAMPLE_ID = "citrobacter_sample"
OUTPUT_DIR = Path(f"output/{SAMPLE_ID}")
ANOMALY_DIR = OUTPUT_DIR / "06_anomaly"

def create_directories():
    """Create anomaly detection directories"""
    dirs = [
        ANOMALY_DIR / "reference",
        ANOMALY_DIR / "alignment", 
        ANOMALY_DIR / "gene_taxonomy"
    ]
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)
    print("Created anomaly directories")


def create_mock_reference_info():
    """Create mock reference selection data"""
    ref_info = {
        "sample_id": SAMPLE_ID,
        "reference": "none",  # No reference available for this test
        "distance": None,
        "status": "no_database"
    }
    
    with open(ANOMALY_DIR / "reference" / "mash_distances.tsv", 'w') as f:
        f.write("# No reference database configured\n")
    
    with open(ANOMALY_DIR / "reference" / "best_reference.txt", 'w') as f:
        f.write("none")
    
    with open(ANOMALY_DIR / "reference" / "reference_info.json", 'w') as f:
        json.dump(ref_info, f, indent=2)
    
    print("Created mock reference info (no database)")


def create_mock_alignment():
    """Create mock alignment data"""
    alignment = {
        "sample_id": SAMPLE_ID,
        "reference": "none",
        "status": "no_reference",
        "aligned_regions": [],
        "unaligned_regions": [],
        "total_aligned_bp": 0,
        "insertions": [],
        "deletions": []
    }
    
    # Create empty coord files
    open(ANOMALY_DIR / "alignment" / "nucmer.delta", 'w').close()
    open(ANOMALY_DIR / "alignment" / "nucmer.coords", 'w').close()
    
    with open(ANOMALY_DIR / "alignment" / "alignment.json", 'w') as f:
        json.dump(alignment, f, indent=2)
    
    print("Created mock alignment data (no reference)")


def create_mock_gene_taxonomy():
    """Create mock gene taxonomy data with some 'foreign' genes for testing"""
    # Simulate finding a few genes from different organisms
    gene_taxonomy = {
        "sample_id": SAMPLE_ID,
        "status": "simulated",
        "host_genus": "Citrobacter",
        "genes": {
            "GENE_001": {"organism": "Citrobacter freundii", "genus": "Citrobacter", "identity": 99.5},
            "GENE_002": {"organism": "Citrobacter freundii", "genus": "Citrobacter", "identity": 98.2},
            "GENE_003": {"organism": "Citrobacter freundii", "genus": "Citrobacter", "identity": 97.8},
            # Simulated foreign genes (for testing visualization)
            "GENE_0456": {"organism": "Salmonella enterica", "genus": "Salmonella", "identity": 95.3},
            "GENE_0457": {"organism": "Salmonella enterica", "genus": "Salmonella", "identity": 94.1},
            "GENE_1234": {"organism": "Escherichia coli", "genus": "Escherichia", "identity": 89.5},
        },
        "foreign_genes": [
            {
                "gene_id": "GENE_0456",
                "organism": "Salmonella enterica",
                "genus": "Salmonella",
                "identity": 95.3
            },
            {
                "gene_id": "GENE_0457",
                "organism": "Salmonella enterica", 
                "genus": "Salmonella",
                "identity": 94.1
            },
            {
                "gene_id": "GENE_1234",
                "organism": "Escherichia coli",
                "genus": "Escherichia",
                "identity": 89.5
            }
        ],
        "taxonomy_distribution": {
            "Citrobacter": 4575,  # Most genes match host
            "Salmonella": 2,
            "Escherichia": 1
        }
    }
    
    # Create empty BLAST results file
    open(ANOMALY_DIR / "gene_taxonomy" / "blast_results.tsv", 'w').close()
    open(ANOMALY_DIR / "gene_taxonomy" / "proteins.faa", 'w').close()
    
    with open(ANOMALY_DIR / "gene_taxonomy" / "gene_taxonomy.json", 'w') as f:
        json.dump(gene_taxonomy, f, indent=2)
    
    print("Created mock gene taxonomy with 3 foreign genes")


def create_anomaly_summary():
    """Create the aggregate anomaly summary"""
    # Load the gene taxonomy
    with open(ANOMALY_DIR / "gene_taxonomy" / "gene_taxonomy.json") as f:
        gene_taxonomy = json.load(f)
    
    with open(ANOMALY_DIR / "reference" / "reference_info.json") as f:
        ref_info = json.load(f)
    
    with open(ANOMALY_DIR / "alignment" / "alignment.json") as f:
        alignment = json.load(f)
    
    summary = {
        "sample_id": SAMPLE_ID,
        
        # Reference info
        "reference": ref_info,
        
        # Alignment summary
        "alignment": {
            "status": alignment.get("status"),
            "total_aligned_bp": alignment.get("total_aligned_bp", 0),
            "num_aligned_regions": len(alignment.get("aligned_regions", [])),
            "insertions": alignment.get("insertions", []),
            "deletions": alignment.get("deletions", [])
        },
        
        # Gene taxonomy
        "gene_taxonomy": {
            "status": gene_taxonomy.get("status"),
            "host_genus": gene_taxonomy.get("host_genus"),
            "foreign_genes": gene_taxonomy.get("foreign_genes", []),
            "num_foreign_genes": len(gene_taxonomy.get("foreign_genes", [])),
            "taxonomy_distribution": gene_taxonomy.get("taxonomy_distribution", {})
        },
        
        # Anomaly flags
        "anomaly_detected": len(gene_taxonomy.get("foreign_genes", [])) > 0,
        "anomaly_types": ["foreign_gene_insertion"] if gene_taxonomy.get("foreign_genes") else []
    }
    
    with open(ANOMALY_DIR / "anomaly_summary.json", 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"Created anomaly summary: {len(gene_taxonomy.get('foreign_genes', []))} foreign genes detected")


def main():
    print("=" * 60)
    print("Generating test anomaly data for visualization")
    print("=" * 60)
    
    create_directories()
    create_mock_reference_info()
    create_mock_alignment()
    create_mock_gene_taxonomy()
    create_anomaly_summary()
    
    print("=" * 60)
    print("Done! Test anomaly data created.")
    print("Run the report generation to see the visualization.")
    print("=" * 60)


if __name__ == "__main__":
    main()

