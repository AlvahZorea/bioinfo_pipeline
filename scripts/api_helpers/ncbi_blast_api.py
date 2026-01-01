#!/usr/bin/env python3
"""
Submit protein sequences to NCBI BLAST API for taxonomic classification.
This is a standalone helper script, NOT part of the Snakemake pipeline.

Usage:
    python ncbi_blast_api.py <input.faa> <output.json> [--max-seqs 100] [--email your@email.com]

The script:
1. Reads protein FASTA
2. Submits batches to NCBI BLAST (blastp vs nr)
3. Parses results to get taxonomic classification per gene
4. Outputs JSON compatible with pipeline format
"""

import argparse
import json
import time
import sys
from pathlib import Path
from collections import defaultdict

try:
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW, NCBIXML
except ImportError:
    print("ERROR: Biopython required. Install with: pip install biopython")
    sys.exit(1)

# Fix SSL certificate issues on Windows
import ssl
import urllib.request

try:
    _create_unverified_https_context = ssl._create_unverified_context
except AttributeError:
    pass
else:
    ssl._create_default_https_context = _create_unverified_https_context


def submit_blast(sequence, program="blastp", database="nr", hitlist_size=5):
    """Submit a single sequence to NCBI BLAST."""
    try:
        result_handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=sequence,
            hitlist_size=hitlist_size,
            expect=1e-10,
            format_type="XML"
        )
        return result_handle
    except Exception as e:
        print(f"  BLAST error: {e}")
        return None


def parse_blast_result(result_handle):
    """Parse BLAST XML result and extract top hit taxonomy."""
    try:
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)
        
        if record.alignments:
            top_hit = record.alignments[0]
            top_hsp = top_hit.hsps[0]
            
            # Parse organism from hit definition
            # Format usually: "gi|xxx|ref|xxx| description [Organism name]"
            title = top_hit.title
            organism = "Unknown"
            if "[" in title and "]" in title:
                organism = title[title.rfind("[")+1:title.rfind("]")]
            
            identity = (top_hsp.identities / top_hsp.align_length) * 100
            
            return {
                "hit_id": top_hit.hit_id,
                "organism": organism,
                "genus": organism.split()[0] if organism else "Unknown",
                "identity": round(identity, 2),
                "evalue": top_hsp.expect,
                "bit_score": top_hsp.bits
            }
    except Exception as e:
        print(f"  Parse error: {e}")
    
    return None


def run_blast_batch(proteins, max_seqs=100, delay=3):
    """Run BLAST on a batch of proteins with rate limiting."""
    
    results = {}
    total = min(len(proteins), max_seqs)
    
    print(f"\nSubmitting {total} proteins to NCBI BLAST...")
    print("(This may take several minutes due to API rate limits)\n")
    
    for i, protein in enumerate(proteins[:max_seqs]):
        seq_id = protein["id"]
        sequence = protein["sequence"]
        
        print(f"[{i+1}/{total}] BLASTing {seq_id}...", end=" ", flush=True)
        
        result_handle = submit_blast(sequence)
        
        if result_handle:
            hit = parse_blast_result(result_handle)
            if hit:
                results[seq_id] = hit
                print(f"-> {hit['organism']} ({hit['identity']:.1f}%)")
            else:
                print("-> No hits")
                results[seq_id] = {"organism": "No hits", "genus": "Unknown", "identity": 0}
        else:
            print("-> Failed")
            results[seq_id] = {"organism": "Error", "genus": "Unknown", "identity": 0}
        
        # Rate limiting - NCBI requires delays between requests
        if i < total - 1:
            time.sleep(delay)
    
    return results


def analyze_taxonomy(results, sample_id="sample"):
    """Analyze BLAST results to identify host and foreign genes."""
    
    # Count genera
    genus_counts = defaultdict(int)
    for seq_id, hit in results.items():
        genus = hit.get("genus", "Unknown")
        if genus and genus != "Unknown" and genus != "No hits":
            genus_counts[genus] += 1
    
    # Determine host (most common genus)
    host_genus = max(genus_counts, key=genus_counts.get) if genus_counts else "Unknown"
    
    # Identify foreign genes
    foreign_genes = []
    for seq_id, hit in results.items():
        genus = hit.get("genus", "Unknown")
        identity = hit.get("identity", 0)
        
        if genus and genus != host_genus and genus != "Unknown" and identity > 70:
            foreign_genes.append({
                "gene_id": seq_id,
                "organism": hit.get("organism", "Unknown"),
                "genus": genus,
                "identity": identity
            })
    
    # Build output
    output = {
        "sample_id": sample_id,
        "status": "completed",
        "tool": "NCBI BLAST API (blastp vs nr)",
        "host_genus": host_genus,
        "genes": results,
        "foreign_genes": foreign_genes,
        "taxonomy_distribution": dict(genus_counts),
        "summary": {
            "total_genes_analyzed": len(results),
            "num_foreign_genes": len(foreign_genes),
            "host_genus": host_genus
        }
    }
    
    return output


def main():
    parser = argparse.ArgumentParser(
        description="Submit proteins to NCBI BLAST for taxonomic classification"
    )
    parser.add_argument("input_fasta", help="Input protein FASTA file")
    parser.add_argument("output_json", help="Output JSON file")
    parser.add_argument("--max-seqs", type=int, default=100, 
                        help="Maximum sequences to BLAST (default: 100)")
    parser.add_argument("--sample-id", default="sample",
                        help="Sample ID for output")
    parser.add_argument("--delay", type=int, default=3,
                        help="Delay between requests in seconds (default: 3)")
    
    args = parser.parse_args()
    
    # Read proteins
    print(f"Reading proteins from {args.input_fasta}...")
    proteins = []
    for record in SeqIO.parse(args.input_fasta, "fasta"):
        proteins.append({
            "id": record.id,
            "sequence": str(record.seq)
        })
    print(f"Found {len(proteins)} proteins")
    
    if len(proteins) == 0:
        print("ERROR: No proteins found in input file")
        sys.exit(1)
    
    # Run BLAST
    results = run_blast_batch(proteins, max_seqs=args.max_seqs, delay=args.delay)
    
    # Analyze results
    output = analyze_taxonomy(results, sample_id=args.sample_id)
    
    # Write output
    with open(args.output_json, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"Results saved to {args.output_json}")
    print(f"{'='*60}")
    print(f"Host genus: {output['host_genus']}")
    print(f"Foreign genes detected: {len(output['foreign_genes'])}")
    
    if output['foreign_genes']:
        print("\nForeign genes:")
        for fg in output['foreign_genes'][:10]:
            print(f"  - {fg['gene_id']}: {fg['organism']} ({fg['identity']:.1f}%)")


if __name__ == "__main__":
    main()

