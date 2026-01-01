#!/usr/bin/env python3
"""
Master script to run all API-based analyses on a sample.
This is a standalone helper script, NOT part of the Snakemake pipeline.

Usage:
    python run_all_analyses.py <sample_dir> [--max-seqs 50] [--skip-blast]

This will:
1. Extract proteins from Bakta annotation
2. Run NCBI BLAST for per-gene taxonomy (optional, slow)
3. Provide instructions for ResFinder/VFDB web submission
4. Update the pipeline JSON files
"""

import argparse
import json
import subprocess
import sys
from pathlib import Path


def run_command(cmd, description):
    """Run a command and return success status."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")
    print(f"Command: {' '.join(cmd)}\n")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False)
        return True
    except subprocess.CalledProcessError as e:
        print(f"ERROR: {e}")
        return False
    except FileNotFoundError:
        print(f"ERROR: Command not found")
        return False


def check_sample_dir(sample_dir):
    """Check that sample directory has required files."""
    required = [
        "02_assembly/contigs.fasta",
        "04_annotation/annotation_stats.json",
    ]
    
    # Check for Bakta output
    bakta_gbff = list(sample_dir.glob("04_annotation/bakta/*.gbff"))
    if not bakta_gbff:
        # Try root directory
        bakta_gbff = list(sample_dir.parent.glob("*.gbff"))
    
    missing = []
    for req in required:
        if not (sample_dir / req).exists():
            missing.append(req)
    
    if not bakta_gbff:
        missing.append("Bakta .gbff file")
    
    return missing, bakta_gbff[0] if bakta_gbff else None


def main():
    parser = argparse.ArgumentParser(
        description="Run all API-based analyses on a sample"
    )
    parser.add_argument("sample_dir", help="Path to sample output directory")
    parser.add_argument("--max-seqs", type=int, default=100,
                        help="Maximum sequences for BLAST (default: 100)")
    parser.add_argument("--skip-blast", action="store_true",
                        help="Skip NCBI BLAST (slow)")
    parser.add_argument("--sample-id", help="Sample ID (default: from directory name)")
    
    args = parser.parse_args()
    
    sample_dir = Path(args.sample_dir)
    sample_id = args.sample_id or sample_dir.name
    
    print("=" * 60)
    print(f"API-Based Analysis Pipeline")
    print(f"Sample: {sample_id}")
    print(f"Directory: {sample_dir}")
    print("=" * 60)
    
    # Check prerequisites
    missing, bakta_gbff = check_sample_dir(sample_dir)
    if missing:
        print(f"\nERROR: Missing required files:")
        for m in missing:
            print(f"  - {m}")
        sys.exit(1)
    
    print(f"\nFound Bakta annotation: {bakta_gbff}")
    
    # Create output directories
    anomaly_dir = sample_dir / "06_anomaly" / "gene_taxonomy"
    anomaly_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Extract proteins
    proteins_fasta = anomaly_dir / "proteins.faa"
    scripts_dir = Path(__file__).parent
    
    success = run_command(
        [sys.executable, str(scripts_dir / "extract_proteins.py"),
         str(bakta_gbff), str(proteins_fasta)],
        "Extract proteins from Bakta annotation"
    )
    
    if not success:
        print("Failed to extract proteins")
        sys.exit(1)
    
    # Step 2: NCBI BLAST for taxonomy (optional)
    if not args.skip_blast:
        taxonomy_json = anomaly_dir / "gene_taxonomy.json"
        
        print("\n" + "!" * 60)
        print("WARNING: NCBI BLAST is SLOW (3+ seconds per protein)")
        print(f"Will analyze up to {args.max_seqs} proteins")
        print("This may take 5-30 minutes depending on --max-seqs")
        print("!" * 60)
        
        response = input("\nProceed with BLAST? [y/N]: ").strip().lower()
        if response == 'y':
            success = run_command(
                [sys.executable, str(scripts_dir / "ncbi_blast_api.py"),
                 str(proteins_fasta), str(taxonomy_json),
                 "--max-seqs", str(args.max_seqs),
                 "--sample-id", sample_id],
                "NCBI BLAST for per-gene taxonomy"
            )
            
            if success:
                # Update anomaly summary
                update_anomaly_summary(sample_dir, taxonomy_json)
        else:
            print("Skipping BLAST")
    else:
        print("\nSkipping NCBI BLAST (--skip-blast)")
    
    # Step 3: Show web submission instructions
    print("\n" + "=" * 60)
    print("WEB-BASED ANALYSES")
    print("=" * 60)
    
    contigs_fasta = sample_dir / "02_assembly" / "contigs.fasta"
    
    print(f"""
For AMR detection, submit to ResFinder:
  1. Go to: https://cge.food.dtu.dk/services/ResFinder/
  2. Upload: {contigs_fasta}
  3. Download results and run:
     python {scripts_dir}/resfinder_api.py {contigs_fasta} \\
         {sample_dir}/05_specialized/amr/amr_summary.json \\
         --resfinder-output <downloaded_results.txt>

For virulence factors, submit to VFDB:
  1. Go to: http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi
  2. Upload: {proteins_fasta}
  3. Or use local BLAST with downloaded VFDB
""")
    
    print("=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"\nTo regenerate the report, run:")
    print(f"  python scripts/generate_report_direct.py")


def update_anomaly_summary(sample_dir, taxonomy_json):
    """Update the anomaly summary with BLAST results."""
    
    anomaly_summary_path = sample_dir / "06_anomaly" / "anomaly_summary.json"
    
    # Load BLAST results
    with open(taxonomy_json) as f:
        taxonomy = json.load(f)
    
    # Load or create anomaly summary
    if anomaly_summary_path.exists():
        with open(anomaly_summary_path) as f:
            summary = json.load(f)
    else:
        summary = {"sample_id": sample_dir.name}
    
    # Update with taxonomy results
    summary["gene_taxonomy"] = {
        "status": taxonomy.get("status"),
        "host_genus": taxonomy.get("host_genus"),
        "foreign_genes": taxonomy.get("foreign_genes", []),
        "num_foreign_genes": len(taxonomy.get("foreign_genes", [])),
        "taxonomy_distribution": taxonomy.get("taxonomy_distribution", {})
    }
    
    # Update anomaly flags
    foreign_genes = taxonomy.get("foreign_genes", [])
    summary["anomaly_detected"] = len(foreign_genes) > 0
    summary["anomaly_types"] = ["foreign_gene_insertion"] if foreign_genes else []
    
    # Save
    with open(anomaly_summary_path, 'w') as f:
        json.dump(summary, f, indent=2)
    
    print(f"\nUpdated anomaly summary: {anomaly_summary_path}")
    print(f"Foreign genes detected: {len(foreign_genes)}")


if __name__ == "__main__":
    main()

