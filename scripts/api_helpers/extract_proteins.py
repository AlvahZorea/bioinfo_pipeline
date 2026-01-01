#!/usr/bin/env python3
"""
Extract protein sequences from Bakta GenBank file for API submissions.
"""

import sys
from pathlib import Path
from Bio import SeqIO


def extract_proteins(gbk_file, output_fasta, max_proteins=None):
    """Extract protein sequences from GenBank file."""
    
    proteins = []
    
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS" and "translation" in feature.qualifiers:
                locus_tag = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                product = feature.qualifiers.get("product", ["hypothetical protein"])[0]
                gene = feature.qualifiers.get("gene", [""])[0]
                translation = feature.qualifiers["translation"][0]
                
                # Get location info
                start = int(feature.location.start)
                end = int(feature.location.end)
                strand = "+" if feature.location.strand == 1 else "-"
                contig = record.id
                
                proteins.append({
                    "id": locus_tag,
                    "gene": gene,
                    "product": product,
                    "sequence": translation,
                    "contig": contig,
                    "start": start,
                    "end": end,
                    "strand": strand
                })
    
    # Limit if specified
    if max_proteins:
        proteins = proteins[:max_proteins]
    
    # Write FASTA
    with open(output_fasta, 'w') as f:
        for p in proteins:
            header = f">{p['id']} {p['gene']} {p['product']} [{p['contig']}:{p['start']}-{p['end']}({p['strand']})]"
            f.write(f"{header}\n{p['sequence']}\n")
    
    print(f"Extracted {len(proteins)} proteins to {output_fasta}")
    return proteins


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python extract_proteins.py <input.gbff> <output.faa> [max_proteins]")
        sys.exit(1)
    
    gbk_file = sys.argv[1]
    output_fasta = sys.argv[2]
    max_proteins = int(sys.argv[3]) if len(sys.argv) > 3 else None
    
    extract_proteins(gbk_file, output_fasta, max_proteins)

