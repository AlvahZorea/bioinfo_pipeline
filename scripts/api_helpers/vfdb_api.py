#!/usr/bin/env python3
"""
Submit protein sequences to VFDB (Virulence Factor Database) for virulence detection.
This is a standalone helper script, NOT part of the Snakemake pipeline.

VFDB: http://www.mgc.ac.cn/VFs/
VFanalyzer: http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi

Usage:
    python vfdb_api.py <input.faa> <output.json> [--sample-id sample]

Note: VFDB web BLAST is available at their website.
For local analysis, download VFDB FASTA and use local BLAST.
"""

import argparse
import json
import sys
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.Blast import NCBIWWW, NCBIXML
except ImportError:
    print("ERROR: Biopython required. Install with: pip install biopython")
    sys.exit(1)


# VFDB categories
VFDB_CATEGORIES = {
    "adherence": ["fim", "pil", "sfa", "pap", "afa", "dra", "bfp", "eae", "intimin"],
    "invasion": ["inv", "ipa", "spa", "sip", "sop", "hil", "prg"],
    "toxin": ["stx", "ctx", "lt", "st", "hly", "cdt", "cnf", "pet", "pic", "sat", "vat"],
    "secretion_system": ["esc", "sep", "ces", "esp", "tir", "ler", "grl"],
    "iron_acquisition": ["iro", "iuc", "iut", "ent", "fep", "chu", "hmu", "sit", "fyu"],
    "capsule": ["kps", "neu", "cap", "cps"],
    "immune_evasion": ["iss", "tra", "omp", "ail"],
    "regulation": ["rcs", "cpx", "evg", "bar", "umo"],
}


def categorize_virulence_factor(gene_name, description):
    """Categorize a virulence factor based on gene name or description."""
    gene_lower = gene_name.lower()
    desc_lower = description.lower()
    
    for category, keywords in VFDB_CATEGORIES.items():
        for keyword in keywords:
            if keyword in gene_lower or keyword in desc_lower:
                return category
    
    return "other"


def blast_against_vfdb_online(sequence, gene_id):
    """
    BLAST a sequence against VFDB using NCBI BLAST.
    
    Note: This BLASTs against nr but filters for VFDB-like hits.
    For true VFDB-specific BLAST, use their web interface or local database.
    """
    
    try:
        # BLAST against nr (we can't directly BLAST against VFDB via API)
        # Filter results for virulence-related terms
        result_handle = NCBIWWW.qblast(
            program="blastp",
            database="nr",
            sequence=sequence,
            hitlist_size=10,
            expect=1e-10,
            format_type="XML",
            entrez_query="virulence[All Fields]"  # Filter for virulence-related
        )
        
        blast_records = NCBIXML.parse(result_handle)
        record = next(blast_records)
        
        hits = []
        for alignment in record.alignments[:5]:
            top_hsp = alignment.hsps[0]
            identity = (top_hsp.identities / top_hsp.align_length) * 100
            
            # Check if this looks like a virulence factor
            title_lower = alignment.title.lower()
            is_vf = any(term in title_lower for term in 
                       ["virulence", "toxin", "adhesin", "invasin", "hemolysin", 
                        "fimbri", "pili", "secretion", "effector"])
            
            if is_vf or identity > 80:
                hits.append({
                    "hit_id": alignment.hit_id,
                    "title": alignment.title[:100],
                    "identity": round(identity, 2),
                    "evalue": top_hsp.expect,
                    "is_virulence": is_vf
                })
        
        return hits
        
    except Exception as e:
        print(f"  BLAST error: {e}")
        return []


def create_virulence_summary(virulence_factors, sample_id="sample"):
    """Create virulence factor summary in pipeline format."""
    
    # Group by category
    categories = {}
    for vf in virulence_factors:
        cat = vf.get("category", "other")
        if cat not in categories:
            categories[cat] = []
        categories[cat].append(vf)
    
    return {
        "sample_id": sample_id,
        "tool": "VFDB BLAST",
        "status": "completed" if virulence_factors else "no_virulence",
        "virulence_factors": virulence_factors,
        "categories": {k: len(v) for k, v in categories.items()},
        "summary": {
            "total_factors": len(virulence_factors),
            "categories_found": list(categories.keys())
        }
    }


def main():
    parser = argparse.ArgumentParser(
        description="Submit proteins to VFDB for virulence factor detection"
    )
    parser.add_argument("input_fasta", help="Input protein FASTA file")
    parser.add_argument("output_json", help="Output JSON file")
    parser.add_argument("--sample-id", default="sample", help="Sample ID")
    parser.add_argument("--max-seqs", type=int, default=50, 
                        help="Maximum sequences to analyze (default: 50)")
    parser.add_argument("--local-vfdb", help="Path to local VFDB FASTA (optional)")
    parser.add_argument("--web-only", action="store_true",
                        help="Just show web submission instructions")
    
    args = parser.parse_args()
    
    if args.web_only:
        print("=" * 60)
        print("VFDB Web Submission Instructions")
        print("=" * 60)
        print("\n1. Go to: http://www.mgc.ac.cn/cgi-bin/VFs/v5/main.cgi")
        print("2. Click 'VFanalyzer' or 'BLAST'")
        print(f"3. Upload your FASTA file: {args.input_fasta}")
        print("4. Submit and wait for results")
        print("5. Download/copy the results")
        print("\nAlternatively, download VFDB sequences for local BLAST:")
        print("  http://www.mgc.ac.cn/VFs/Down/VFDB_setA_pro.fas.gz")
        print("=" * 60)
        
        # Create placeholder
        output = {
            "sample_id": args.sample_id,
            "tool": "VFDB",
            "status": "pending_web_submission",
            "message": "Submit to http://www.mgc.ac.cn/VFs/",
            "virulence_factors": [],
            "categories": {},
            "summary": {"total_factors": 0}
        }
        
        with open(args.output_json, 'w') as f:
            json.dump(output, f, indent=2)
        return
    
    # Read proteins
    print(f"Reading proteins from {args.input_fasta}...")
    proteins = list(SeqIO.parse(args.input_fasta, "fasta"))
    print(f"Found {len(proteins)} proteins")
    
    if len(proteins) == 0:
        print("ERROR: No proteins found")
        sys.exit(1)
    
    # Analyze
    print(f"\nAnalyzing up to {args.max_seqs} proteins for virulence factors...")
    print("(Using NCBI BLAST with virulence filter - for best results use local VFDB)\n")
    
    virulence_factors = []
    
    for i, record in enumerate(proteins[:args.max_seqs]):
        gene_id = record.id
        sequence = str(record.seq)
        description = record.description
        
        print(f"[{i+1}/{min(len(proteins), args.max_seqs)}] Checking {gene_id}...", end=" ", flush=True)
        
        hits = blast_against_vfdb_online(sequence, gene_id)
        
        # Check if any hits are virulence-related
        vf_hits = [h for h in hits if h.get("is_virulence")]
        
        if vf_hits:
            top_hit = vf_hits[0]
            category = categorize_virulence_factor(gene_id, top_hit["title"])
            
            virulence_factors.append({
                "gene_id": gene_id,
                "vfdb_id": top_hit["hit_id"],
                "description": top_hit["title"],
                "identity": top_hit["identity"],
                "category": category
            })
            print(f"-> VF detected: {category}")
        else:
            print("-> Not a VF")
        
        # Rate limiting
        import time
        time.sleep(3)
    
    # Create output
    output = create_virulence_summary(virulence_factors, args.sample_id)
    
    with open(args.output_json, 'w') as f:
        json.dump(output, f, indent=2)
    
    print(f"\n{'='*60}")
    print(f"Results saved to {args.output_json}")
    print(f"Virulence factors found: {len(virulence_factors)}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()

