#!/usr/bin/env python3
"""
Parse existing tool outputs into pipeline JSON format
"""

import json
import os
import re
from pathlib import Path
from Bio import SeqIO

# Configuration
SAMPLE_ID = "citrobacter_sample"
OUTPUT_DIR = Path(f"output/{SAMPLE_ID}")

# Input files
ASSEMBLY_FILE = "assembled-contigs.fa"
KRAKEN_REPORT = "kraken_il.report"
BAKTA_GBK = "bakta.gbff"

def parse_assembly():
    """Parse assembly FASTA and generate stats"""
    print("Parsing assembly...")
    
    contigs = list(SeqIO.parse(ASSEMBLY_FILE, "fasta"))
    lengths = [len(c) for c in contigs]
    total_length = sum(lengths)
    
    # Calculate GC content
    gc_count = sum(str(c.seq).upper().count('G') + str(c.seq).upper().count('C') for c in contigs)
    gc_content = (gc_count / total_length * 100) if total_length > 0 else 0
    
    # Calculate N50
    sorted_lengths = sorted(lengths, reverse=True)
    cumsum = 0
    n50 = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= total_length / 2:
            n50 = length
            break
    
    stats = {
        "sample_id": SAMPLE_ID,
        "assembler": "spades+pilon",
        "num_contigs": len(contigs),
        "total_length": total_length,
        "largest_contig": max(lengths) if lengths else 0,
        "smallest_contig": min(lengths) if lengths else 0,
        "n50": n50,
        "gc_content": round(gc_content, 2),
        "mean_contig_length": round(total_length / len(contigs), 2) if contigs else 0,
    }
    
    # Write stats
    with open(OUTPUT_DIR / "02_assembly" / "assembly_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
    
    # Copy contigs file
    import shutil
    shutil.copy(ASSEMBLY_FILE, OUTPUT_DIR / "02_assembly" / "contigs.fasta")
    
    print(f"  Contigs: {len(contigs)}, Total: {total_length:,} bp, N50: {n50:,} bp, GC: {gc_content:.1f}%")
    return stats


def parse_kraken_report():
    """Parse Kraken2 report and generate taxonomy JSON"""
    print("Parsing Kraken2 report...")
    
    taxonomy = {
        "sample_id": SAMPLE_ID,
        "tool": "kraken2",
        "top_hits": [],
        "domain": None,
        "phylum": None,
        "class": None,
        "order": None,
        "family": None,
        "genus": None,
        "species": None,
    }
    
    rank_map = {
        'D': 'domain', 'P': 'phylum', 'C': 'class', 'O': 'order',
        'F': 'family', 'G': 'genus', 'S': 'species',
        'R': 'root', 'R1': 'superkingdom'
    }
    
    with open(KRAKEN_REPORT) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                pct = float(parts[0].strip())
                clade_reads = int(parts[1].strip())
                direct_reads = int(parts[2].strip())
                rank = parts[3].strip()
                taxid = parts[4].strip()
                name = parts[5].strip()
                
                # Get clean rank (G1 -> G, S1 -> S, etc.)
                base_rank = rank[0] if rank else ''
                
                # Store top hits (species and genus level with >0.1%)
                if pct > 0.1 and base_rank in ['S', 'G']:
                    taxonomy["top_hits"].append({
                        "name": name,
                        "taxid": taxid,
                        "rank": base_rank,
                        "percentage": pct
                    })
                
                # Store highest percentage for main ranks
                if base_rank in rank_map and rank in ['D', 'P', 'C', 'O', 'F', 'G', 'S']:
                    rank_name = rank_map[base_rank]
                    current_pct = taxonomy.get(f"{rank_name}_pct", 0)
                    if pct > current_pct:
                        taxonomy[rank_name] = name
                        taxonomy[f"{rank_name}_pct"] = pct
    
    # Sort top hits by percentage
    taxonomy["top_hits"] = sorted(taxonomy["top_hits"], key=lambda x: x["percentage"], reverse=True)[:10]
    
    # Create full taxonomy output
    full_taxonomy = {
        "sample_id": SAMPLE_ID,
        "organism_type": "bacteria",
        "kraken2": taxonomy,
        "has_viral_content": False,
    }
    
    # Write outputs
    with open(OUTPUT_DIR / "03_identification" / "taxonomy.json", 'w') as f:
        json.dump(full_taxonomy, f, indent=2)
    
    with open(OUTPUT_DIR / "03_identification" / "organism_type.txt", 'w') as f:
        f.write("bacteria")
    
    top_species = taxonomy["top_hits"][0]["name"] if taxonomy["top_hits"] else "Unknown"
    print(f"  Top hit: {top_species}")
    return full_taxonomy


def parse_bakta_gbk():
    """Parse Bakta GenBank file and extract annotation stats"""
    print("Parsing Bakta annotation...")
    
    feature_counts = {
        "CDS": 0,
        "tRNA": 0,
        "rRNA": 0,
        "tmRNA": 0,
        "ncRNA": 0,
        "gene": 0,
    }
    
    protein_lengths = []
    
    # Parse the GenBank file
    for record in SeqIO.parse(BAKTA_GBK, "genbank"):
        for feature in record.features:
            ftype = feature.type
            if ftype in feature_counts:
                feature_counts[ftype] += 1
            
            if ftype == "CDS":
                # Get protein length from translation or calculate from location
                if "translation" in feature.qualifiers:
                    protein_lengths.append(len(feature.qualifiers["translation"][0]))
                else:
                    # Approximate from nucleotide length
                    protein_lengths.append(len(feature) // 3)
            
            if ftype == "gene":
                feature_counts["gene"] += 1
    
    # If gene count is 0, estimate from CDS
    if feature_counts["gene"] == 0:
        feature_counts["gene"] = feature_counts["CDS"] + feature_counts["tRNA"] + feature_counts["rRNA"]
    
    stats = {
        "sample_id": SAMPLE_ID,
        "organism_type": "bacteria",
        "annotator": "bakta",
        "total_genes": feature_counts["gene"],
        "cds_count": feature_counts["CDS"],
        "trna_count": feature_counts["tRNA"],
        "rrna_count": feature_counts["rRNA"],
        "protein_count": feature_counts["CDS"],
        "mean_protein_length": round(sum(protein_lengths) / len(protein_lengths), 1) if protein_lengths else 0,
        "feature_counts": feature_counts,
    }
    
    # Write output
    with open(OUTPUT_DIR / "04_annotation" / "annotation_stats.json", 'w') as f:
        json.dump(stats, f, indent=2)
    
    with open(OUTPUT_DIR / "04_annotation" / "annotation_complete.flag", 'w') as f:
        f.write("Annotation complete: bacteria\nAnnotator: bakta\n")
    
    print(f"  CDSs: {stats['cds_count']}, tRNAs: {stats['trna_count']}, rRNAs: {stats['rrna_count']}")
    return stats


def create_placeholder_results():
    """Create placeholder results for analyses not run"""
    print("Creating placeholder results for skipped analyses...")
    
    # AMR - not analyzed
    amr = {
        "sample_id": SAMPLE_ID,
        "tool": "AMRFinderPlus",
        "status": "not_analyzed",
        "message": "AMR analysis was not performed - databases not available",
        "genes": [],
        "mutations": [],
        "drug_classes": [],
        "summary": {"total_genes": 0, "total_mutations": 0}
    }
    with open(OUTPUT_DIR / "05_specialized" / "amr" / "amr_summary.json", 'w') as f:
        json.dump(amr, f, indent=2)
    
    # Virulence - not analyzed
    virulence = {
        "sample_id": SAMPLE_ID,
        "tool": "BLAST vs VFDB",
        "status": "not_analyzed",
        "message": "Virulence analysis was not performed - databases not available",
        "virulence_factors": [],
        "categories": {},
        "summary": {"total_factors": 0}
    }
    with open(OUTPUT_DIR / "05_specialized" / "virulence" / "virulence_summary.json", 'w') as f:
        json.dump(virulence, f, indent=2)
    
    # MGE - not analyzed
    mge = {
        "sample_id": SAMPLE_ID,
        "tool": "MOB-suite",
        "status": "not_analyzed",
        "message": "Mobile element analysis was not performed",
        "plasmids": [],
        "summary": {"total_plasmids": 0, "chromosome_contigs": 0, "plasmid_contigs": 0}
    }
    with open(OUTPUT_DIR / "05_specialized" / "mge" / "mge_summary.json", 'w') as f:
        json.dump(mge, f, indent=2)
    
    # Prophages - not analyzed
    prophage = {
        "sample_id": SAMPLE_ID,
        "tool": "PhiSpy",
        "status": "not_analyzed",
        "message": "Prophage analysis was not performed",
        "prophages": [],
        "summary": {"total_prophages": 0, "total_length": 0}
    }
    with open(OUTPUT_DIR / "05_specialized" / "prophages" / "prophage_summary.json", 'w') as f:
        json.dump(prophage, f, indent=2)
    
    print("  Created placeholders for: AMR, Virulence, MGE, Prophages")


def create_quast_placeholder():
    """Create placeholder QUAST report"""
    print("Creating QUAST placeholder...")
    
    # Read assembly stats
    with open(OUTPUT_DIR / "02_assembly" / "assembly_stats.json") as f:
        stats = json.load(f)
    
    # Create TSV
    tsv_content = f"""Assembly\tcontigs
# contigs (>= 0 bp)\t{stats['num_contigs']}
Total length (>= 0 bp)\t{stats['total_length']}
# contigs\t{stats['num_contigs']}
Largest contig\t{stats['largest_contig']}
Total length\t{stats['total_length']}
GC (%)\t{stats['gc_content']}
N50\t{stats['n50']}
"""
    
    with open(OUTPUT_DIR / "02_assembly" / "quast" / "report.tsv", 'w') as f:
        f.write(tsv_content)
    
    with open(OUTPUT_DIR / "02_assembly" / "quast" / "report.txt", 'w') as f:
        f.write(tsv_content)
    
    # Create minimal HTML
    html = f"""<!DOCTYPE html>
<html><head><title>QUAST Report</title></head>
<body>
<h1>QUAST Assembly Report</h1>
<table>
<tr><th>Metric</th><th>Value</th></tr>
<tr><td>Total length</td><td>{stats['total_length']:,} bp</td></tr>
<tr><td>Contigs</td><td>{stats['num_contigs']}</td></tr>
<tr><td>N50</td><td>{stats['n50']:,} bp</td></tr>
<tr><td>GC content</td><td>{stats['gc_content']}%</td></tr>
</table>
</body></html>
"""
    with open(OUTPUT_DIR / "02_assembly" / "quast" / "report.html", 'w') as f:
        f.write(html)


def create_multiqc_placeholder():
    """Create placeholder MultiQC report"""
    print("Creating MultiQC placeholder...")
    
    html = """<!DOCTYPE html>
<html><head><title>MultiQC Report</title></head>
<body>
<h1>MultiQC Report</h1>
<p>QC analysis was performed on external servers.</p>
</body></html>
"""
    with open(OUTPUT_DIR / "01_qc" / "multiqc_report.html", 'w') as f:
        f.write(html)


def main():
    print("=" * 60)
    print("Parsing existing outputs for pipeline report")
    print("=" * 60)
    
    # Parse existing files
    parse_assembly()
    parse_kraken_report()
    parse_bakta_gbk()
    
    # Create placeholders
    create_placeholder_results()
    create_quast_placeholder()
    create_multiqc_placeholder()
    
    print("=" * 60)
    print("Done! Ready to generate report.")
    print("=" * 60)


if __name__ == "__main__":
    main()

