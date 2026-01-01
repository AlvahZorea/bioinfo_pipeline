"""
============================================================================
Module 6: Anomaly Detection
============================================================================
Detect genomic editing events:
- Foreign gene insertions (via per-gene BLAST taxonomy)
- Deletions/insertions (via reference alignment)
- Integrated plasmids
- Mosaic genomes

Does NOT include: GC/codon compositional analysis
"""

import os

# ============================================================================
# Reference Selection (uses Mash)
# ============================================================================

rule mash_sketch:
    """Create Mash sketch of assembly for reference matching"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta"
    output:
        sketch = f"{OUTPUT_DIR}/06_anomaly/reference/{SAMPLE_ID}.msh"
    params:
        outdir = f"{OUTPUT_DIR}/06_anomaly/reference"
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/mash_sketch.log"
    conda:
        "../envs/anomaly.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        mash sketch -o {params.outdir}/{SAMPLE_ID} -p {threads} {input.assembly} 2>&1 | tee {log}
        """


rule select_reference:
    """Find closest reference genome using Mash"""
    input:
        sketch = rules.mash_sketch.output.sketch
    output:
        distances = f"{OUTPUT_DIR}/06_anomaly/reference/mash_distances.tsv",
        best_ref = f"{OUTPUT_DIR}/06_anomaly/reference/best_reference.txt",
        ref_info = f"{OUTPUT_DIR}/06_anomaly/reference/reference_info.json"
    params:
        refseq_sketches = lambda w: config.get("databases", {}).get("refseq_mash", ""),
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/mash_select_ref.log"
    conda:
        "../envs/anomaly.yaml"
    run:
        import json
        import subprocess
        
        ref_info = {
            "sample_id": SAMPLE_ID,
            "reference": None,
            "distance": None,
            "status": "no_database"
        }
        
        refseq_path = params.refseq_sketches
        
        if refseq_path and os.path.exists(refseq_path):
            # Run mash dist
            cmd = f"mash dist -p {threads} {input.sketch} {refseq_path} > {output.distances} 2>> {log}"
            subprocess.run(cmd, shell=True)
            
            # Find best match
            best_dist = float('inf')
            best_ref = "none"
            
            if os.path.exists(output.distances):
                with open(output.distances) as f:
                    for line in f:
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            ref = parts[1]
                            dist = float(parts[2])
                            if dist < best_dist:
                                best_dist = dist
                                best_ref = ref
            
            if best_ref != "none":
                ref_info["reference"] = best_ref
                ref_info["distance"] = best_dist
                ref_info["status"] = "found"
        else:
            # No database - create empty outputs
            with open(output.distances, 'w') as f:
                f.write("# No reference database configured\n")
        
        with open(output.best_ref, 'w') as f:
            f.write(ref_info.get("reference", "none") or "none")
        
        with open(output.ref_info, 'w') as f:
            json.dump(ref_info, f, indent=2)


# ============================================================================
# Reference Alignment (MUMmer)
# ============================================================================

rule nucmer_align:
    """Align assembly to reference genome using nucmer"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        best_ref = f"{OUTPUT_DIR}/06_anomaly/reference/best_reference.txt"
    output:
        delta = f"{OUTPUT_DIR}/06_anomaly/alignment/nucmer.delta",
        coords = f"{OUTPUT_DIR}/06_anomaly/alignment/nucmer.coords",
        alignment_json = f"{OUTPUT_DIR}/06_anomaly/alignment/alignment.json"
    params:
        outdir = f"{OUTPUT_DIR}/06_anomaly/alignment",
        prefix = f"{OUTPUT_DIR}/06_anomaly/alignment/nucmer"
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/nucmer.log"
    conda:
        "../envs/anomaly.yaml"
    run:
        import json
        import subprocess
        import os
        
        os.makedirs(params.outdir, exist_ok=True)
        
        # Read reference path
        with open(input.best_ref) as f:
            ref_path = f.read().strip()
        
        alignment = {
            "sample_id": SAMPLE_ID,
            "reference": ref_path,
            "status": "no_reference",
            "aligned_regions": [],
            "unaligned_regions": [],
            "total_aligned_bp": 0,
            "insertions": [],
            "deletions": []
        }
        
        if ref_path != "none" and os.path.exists(ref_path):
            # Run nucmer
            cmd = f"nucmer --prefix={params.prefix} -t {threads} {ref_path} {input.assembly} 2>&1 | tee {log}"
            subprocess.run(cmd, shell=True)
            
            # Generate coords
            if os.path.exists(output.delta):
                cmd = f"show-coords -rcl {output.delta} > {output.coords} 2>> {log}"
                subprocess.run(cmd, shell=True)
                alignment["status"] = "aligned"
                
                # Parse coords to get aligned regions
                if os.path.exists(output.coords):
                    with open(output.coords) as f:
                        for line in f:
                            if line.startswith('=') or line.startswith('[') or '|' in line:
                                continue
                            parts = line.strip().split()
                            if len(parts) >= 12:
                                try:
                                    alignment["aligned_regions"].append({
                                        "ref_start": int(parts[0]),
                                        "ref_end": int(parts[1]),
                                        "qry_start": int(parts[3]),
                                        "qry_end": int(parts[4]),
                                        "identity": float(parts[9]),
                                        "contig": parts[11]
                                    })
                                    alignment["total_aligned_bp"] += abs(int(parts[1]) - int(parts[0]))
                                except (ValueError, IndexError):
                                    pass
        else:
            # Create empty files
            open(output.delta, 'w').close()
            open(output.coords, 'w').close()
        
        with open(output.alignment_json, 'w') as f:
            json.dump(alignment, f, indent=2)


# ============================================================================
# Per-Gene BLAST Taxonomy (uses existing Bakta proteins)
# ============================================================================

rule extract_proteins:
    """Extract protein sequences from Bakta annotation (if not already available)"""
    input:
        annotation_flag = f"{OUTPUT_DIR}/04_annotation/annotation_complete.flag"
    output:
        proteins = f"{OUTPUT_DIR}/06_anomaly/gene_taxonomy/proteins.faa"
    params:
        bakta_proteins = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.faa",
        gbk_file = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gbff"
    run:
        import os
        import shutil
        from Bio import SeqIO
        
        os.makedirs(os.path.dirname(output.proteins), exist_ok=True)
        
        # Try to use existing Bakta proteins
        if os.path.exists(params.bakta_proteins):
            shutil.copy(params.bakta_proteins, output.proteins)
        elif os.path.exists(params.gbk_file):
            # Extract from GenBank file
            with open(output.proteins, 'w') as out_f:
                for record in SeqIO.parse(params.gbk_file, "genbank"):
                    for feature in record.features:
                        if feature.type == "CDS" and "translation" in feature.qualifiers:
                            locus = feature.qualifiers.get("locus_tag", ["unknown"])[0]
                            product = feature.qualifiers.get("product", ["hypothetical protein"])[0]
                            seq = feature.qualifiers["translation"][0]
                            out_f.write(f">{locus} {product}\n{seq}\n")
        else:
            # Create empty file
            open(output.proteins, 'w').close()


rule blast_genes:
    """BLAST genes against nr/nt to determine taxonomic origin of each gene"""
    input:
        proteins = rules.extract_proteins.output.proteins
    output:
        blast_out = f"{OUTPUT_DIR}/06_anomaly/gene_taxonomy/blast_results.tsv",
        gene_taxonomy = f"{OUTPUT_DIR}/06_anomaly/gene_taxonomy/gene_taxonomy.json"
    params:
        blast_db = lambda w: config.get("databases", {}).get("nr", "") or config.get("databases", {}).get("refseq_protein", ""),
        evalue = 1e-10,
        max_targets = 3
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/blast_genes.log"
    conda:
        "../envs/anomaly.yaml"
    run:
        import json
        import subprocess
        import os
        from collections import defaultdict
        
        gene_taxonomy = {
            "sample_id": SAMPLE_ID,
            "status": "not_run",
            "host_genus": None,
            "genes": {},
            "foreign_genes": [],
            "taxonomy_distribution": {}
        }
        
        blast_db = params.blast_db
        
        # Check if we have proteins and database
        if os.path.getsize(input.proteins) > 0 and blast_db and os.path.exists(blast_db + ".phr"):
            # Run BLASTP
            cmd = f"""blastp \
                -query {input.proteins} \
                -db {blast_db} \
                -out {output.blast_out} \
                -outfmt "6 qseqid sseqid pident length evalue staxids sscinames" \
                -evalue {params.evalue} \
                -max_target_seqs {params.max_targets} \
                -num_threads {threads} 2>&1 | tee {log}"""
            subprocess.run(cmd, shell=True)
            
            gene_taxonomy["status"] = "completed"
            
            # Parse results
            genus_counts = defaultdict(int)
            
            with open(output.blast_out) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 7:
                        gene_id = parts[0]
                        pident = float(parts[2])
                        organism = parts[6] if len(parts) > 6 else "Unknown"
                        genus = organism.split()[0] if organism else "Unknown"
                        
                        if gene_id not in gene_taxonomy["genes"]:
                            gene_taxonomy["genes"][gene_id] = {
                                "organism": organism,
                                "genus": genus,
                                "identity": pident
                            }
                            genus_counts[genus] += 1
            
            # Determine host (most common genus)
            if genus_counts:
                host = max(genus_counts, key=genus_counts.get)
                gene_taxonomy["host_genus"] = host
                gene_taxonomy["taxonomy_distribution"] = dict(genus_counts)
                
                # Flag foreign genes
                for gene_id, data in gene_taxonomy["genes"].items():
                    if data["genus"] != host and data["identity"] > 70:
                        gene_taxonomy["foreign_genes"].append({
                            "gene_id": gene_id,
                            "organism": data["organism"],
                            "genus": data["genus"],
                            "identity": data["identity"]
                        })
        else:
            # No BLAST db available
            gene_taxonomy["status"] = "no_database"
            open(output.blast_out, 'w').close()
        
        with open(output.gene_taxonomy, 'w') as f:
            json.dump(gene_taxonomy, f, indent=2)


# ============================================================================
# Aggregate Anomaly Results
# ============================================================================

rule aggregate_anomalies:
    """Combine all anomaly detection results"""
    input:
        ref_info = f"{OUTPUT_DIR}/06_anomaly/reference/reference_info.json",
        alignment = f"{OUTPUT_DIR}/06_anomaly/alignment/alignment.json",
        gene_taxonomy = f"{OUTPUT_DIR}/06_anomaly/gene_taxonomy/gene_taxonomy.json"
    output:
        anomaly_summary = f"{OUTPUT_DIR}/06_anomaly/anomaly_summary.json"
    run:
        import json
        
        # Load all inputs
        with open(input.ref_info) as f:
            ref_info = json.load(f)
        
        with open(input.alignment) as f:
            alignment = json.load(f)
        
        with open(input.gene_taxonomy) as f:
            gene_taxonomy = json.load(f)
        
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
            "anomaly_detected": False,
            "anomaly_types": []
        }
        
        # Check for anomalies
        if summary["gene_taxonomy"]["num_foreign_genes"] > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("foreign_gene_insertion")
        
        if len(summary["alignment"]["insertions"]) > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("structural_insertion")
        
        if len(summary["alignment"]["deletions"]) > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("structural_deletion")
        
        with open(output.anomaly_summary, 'w') as f:
            json.dump(summary, f, indent=2)
