#!/usr/bin/env python3
"""
Submit sequences to ResFinder web service for AMR gene detection.
This is a standalone helper script, NOT part of the Snakemake pipeline.

ResFinder is provided by the Center for Genomic Epidemiology (CGE).
Web service: https://cge.food.dtu.dk/services/ResFinder/

Usage:
    python resfinder_api.py <input.fasta> <output.json> [--species "Escherichia coli"]

Note: This uses the CGE web API which has usage limits.
For heavy usage, consider installing ResFinder locally.
"""

import argparse
import json
import requests
import time
import sys
from pathlib import Path

# CGE ResFinder API endpoint
RESFINDER_URL = "https://cge.food.dtu.dk/services/ResFinder/api/"


def read_fasta(fasta_file):
    """Read FASTA file and return as single string."""
    sequences = []
    current_seq = []
    current_header = None
    
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_header and current_seq:
                    sequences.append((current_header, "".join(current_seq)))
                current_header = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_header and current_seq:
            sequences.append((current_header, "".join(current_seq)))
    
    return sequences


def submit_to_resfinder(fasta_content, species=None):
    """
    Submit sequence to ResFinder API.
    
    Note: The CGE API may require registration or have rate limits.
    This is a simplified implementation - check CGE docs for current API specs.
    """
    
    # For now, we'll use a simpler approach with their web form
    # The actual API may require authentication
    
    print("Note: ResFinder web API may require manual submission at:")
    print("  https://cge.food.dtu.dk/services/ResFinder/")
    print("\nAlternatively, you can use the command-line version:")
    print("  pip install resfinder")
    print("  python -m resfinder -ifa input.fasta -o output_dir")
    
    return None


def parse_resfinder_results(result_file):
    """Parse ResFinder output (if run locally or downloaded from web)."""
    
    amr_genes = []
    
    # ResFinder outputs a tab-separated file
    if Path(result_file).exists():
        with open(result_file) as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 8:
                    amr_genes.append({
                        "gene_symbol": parts[0],
                        "identity": float(parts[1]),
                        "coverage": float(parts[2]),
                        "contig": parts[5],
                        "position": parts[6],
                        "accession": parts[7],
                        "drug_class": parts[8] if len(parts) > 8 else "Unknown"
                    })
    
    return amr_genes


def create_amr_summary(amr_genes, sample_id="sample"):
    """Create AMR summary in pipeline format."""
    
    drug_classes = list(set(g.get("drug_class", "Unknown") for g in amr_genes))
    
    return {
        "sample_id": sample_id,
        "tool": "ResFinder",
        "status": "completed" if amr_genes else "no_resistance",
        "genes": amr_genes,
        "mutations": [],  # ResFinder focuses on acquired genes
        "drug_classes": drug_classes,
        "summary": {
            "total_genes": len(amr_genes),
            "total_mutations": 0,
            "drug_classes_affected": len(drug_classes)
        }
    }


def main():
    parser = argparse.ArgumentParser(
        description="Submit sequences to ResFinder for AMR detection"
    )
    parser.add_argument("input_fasta", help="Input nucleotide FASTA file")
    parser.add_argument("output_json", help="Output JSON file")
    parser.add_argument("--species", help="Species name (optional)")
    parser.add_argument("--sample-id", default="sample", help="Sample ID")
    parser.add_argument("--resfinder-output", help="Path to ResFinder results file (if already run)")
    
    args = parser.parse_args()
    
    if args.resfinder_output and Path(args.resfinder_output).exists():
        # Parse existing ResFinder output
        print(f"Parsing ResFinder results from {args.resfinder_output}...")
        amr_genes = parse_resfinder_results(args.resfinder_output)
        
        output = create_amr_summary(amr_genes, args.sample_id)
        
        with open(args.output_json, 'w') as f:
            json.dump(output, f, indent=2)
        
        print(f"Found {len(amr_genes)} AMR genes")
        print(f"Results saved to {args.output_json}")
    else:
        # Show instructions for web submission
        print("=" * 60)
        print("ResFinder Web Submission Instructions")
        print("=" * 60)
        print(f"\n1. Go to: https://cge.food.dtu.dk/services/ResFinder/")
        print(f"2. Upload your FASTA file: {args.input_fasta}")
        if args.species:
            print(f"3. Select species: {args.species}")
        print("4. Submit and wait for results")
        print("5. Download the results file")
        print(f"6. Run: python {sys.argv[0]} {args.input_fasta} {args.output_json} --resfinder-output <results.txt>")
        print("\n" + "=" * 60)
        
        # Create placeholder output
        output = {
            "sample_id": args.sample_id,
            "tool": "ResFinder",
            "status": "pending_web_submission",
            "message": "Submit to https://cge.food.dtu.dk/services/ResFinder/",
            "genes": [],
            "mutations": [],
            "drug_classes": [],
            "summary": {"total_genes": 0, "total_mutations": 0}
        }
        
        with open(args.output_json, 'w') as f:
            json.dump(output, f, indent=2)


if __name__ == "__main__":
    main()

