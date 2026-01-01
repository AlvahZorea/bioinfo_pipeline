"""
============================================================================
Module 6: Anomaly Detection
============================================================================
Detect genomic editing events:
- Foreign gene insertions (via per-gene BLAST taxonomy)
- Deletions/insertions (via reference alignment)
- Integrated plasmids
- Mosaic genomes
- Gene rearrangements (synteny analysis)
- Foreign regulatory elements (intergenic scanning)
- Novel/synthetic genes (no BLAST hits)

Does NOT include: GC/codon compositional analysis
"""

import os
from collections import defaultdict

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


rule download_reference_ncbi:
    """Auto-download reference genome from NCBI based on Kraken2 taxonomy"""
    input:
        taxonomy = f"{OUTPUT_DIR}/03_identification/taxonomy.json",
        ref_info = f"{OUTPUT_DIR}/06_anomaly/reference/reference_info.json"
    output:
        ref_fasta = f"{OUTPUT_DIR}/06_anomaly/reference/ncbi_reference.fasta",
        ref_gff = f"{OUTPUT_DIR}/06_anomaly/reference/ncbi_reference.gff",
        download_info = f"{OUTPUT_DIR}/06_anomaly/reference/download_info.json"
    params:
        outdir = f"{OUTPUT_DIR}/06_anomaly/reference",
        cache_dir = lambda w: config.get("databases", {}).get("ncbi_cache", "databases/ncbi_cache")
    log:
        f"{OUTPUT_DIR}/00_logs/ncbi_download.log"
    run:
        import json
        import subprocess
        import urllib.request
        import os
        import gzip
        import shutil
        
        os.makedirs(params.outdir, exist_ok=True)
        os.makedirs(params.cache_dir, exist_ok=True)
        
        download_info = {
            "sample_id": SAMPLE_ID,
            "status": "not_downloaded",
            "species": None,
            "accession": None,
            "source": None
        }
        
        # Check if we already have a local reference
        with open(input.ref_info) as f:
            ref_info = json.load(f)
        
        if ref_info.get("status") == "found" and ref_info.get("reference"):
            # Already have local reference, just copy/link
            download_info["status"] = "using_local"
            download_info["source"] = ref_info["reference"]
            
            # Create symlinks or copy
            if os.path.exists(ref_info["reference"]):
                shutil.copy(ref_info["reference"], output.ref_fasta)
                # Try to find associated GFF
                gff_path = ref_info["reference"].replace(".fasta", ".gff").replace(".fna", ".gff")
                if os.path.exists(gff_path):
                    shutil.copy(gff_path, output.ref_gff)
                else:
                    open(output.ref_gff, 'w').close()
        else:
            # Need to download from NCBI
            with open(input.taxonomy) as f:
                taxonomy = json.load(f)
            
            # Get species from Kraken2 results
            kraken = taxonomy.get("kraken2", {})
            species = kraken.get("species") or ""
            
            if not species and kraken.get("top_hits"):
                for hit in kraken["top_hits"]:
                    if hit.get("rank") == "S":
                        species = hit.get("name", "")
                        break
            
            download_info["species"] = species
            
            if species:
                # Try to download using NCBI datasets
                species_query = species.replace(" ", "_")
                cache_file = os.path.join(params.cache_dir, f"{species_query}.fasta")
                cache_gff = os.path.join(params.cache_dir, f"{species_query}.gff")
                
                # Check cache first
                if os.path.exists(cache_file):
                    shutil.copy(cache_file, output.ref_fasta)
                    if os.path.exists(cache_gff):
                        shutil.copy(cache_gff, output.ref_gff)
                    else:
                        open(output.ref_gff, 'w').close()
                    download_info["status"] = "from_cache"
                    download_info["source"] = cache_file
                else:
                    # Try NCBI Datasets API
                    try:
                        # Use datasets command if available
                        cmd = f"""datasets download genome taxon "{species}" \
                            --reference --include genome,gff3 \
                            --filename {params.outdir}/ncbi_download.zip 2>&1 | tee {log}"""
                        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
                        
                        if result.returncode == 0 and os.path.exists(f"{params.outdir}/ncbi_download.zip"):
                            # Extract
                            subprocess.run(f"unzip -o {params.outdir}/ncbi_download.zip -d {params.outdir}/ncbi_extract", shell=True)
                            
                            # Find and copy files
                            for root, dirs, files in os.walk(f"{params.outdir}/ncbi_extract"):
                                for f in files:
                                    if f.endswith(".fna"):
                                        shutil.copy(os.path.join(root, f), output.ref_fasta)
                                        shutil.copy(os.path.join(root, f), cache_file)
                                    elif f.endswith(".gff"):
                                        shutil.copy(os.path.join(root, f), output.ref_gff)
                                        shutil.copy(os.path.join(root, f), cache_gff)
                            
                            download_info["status"] = "downloaded"
                            download_info["source"] = "NCBI Datasets"
                        else:
                            # Fallback: create empty files
                            open(output.ref_fasta, 'w').close()
                            open(output.ref_gff, 'w').close()
                            download_info["status"] = "download_failed"
                            download_info["error"] = result.stderr[:500] if result.stderr else "Unknown error"
                    except Exception as e:
                        open(output.ref_fasta, 'w').close()
                        open(output.ref_gff, 'w').close()
                        download_info["status"] = "error"
                        download_info["error"] = str(e)[:500]
            else:
                # No species identified
                open(output.ref_fasta, 'w').close()
                open(output.ref_gff, 'w').close()
                download_info["status"] = "no_species_identified"
        
        with open(output.download_info, 'w') as f:
            json.dump(download_info, f, indent=2)


# ============================================================================
# Reference Alignment (MUMmer)
# ============================================================================

rule nucmer_align:
    """Align assembly to reference genome using nucmer"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        best_ref = f"{OUTPUT_DIR}/06_anomaly/reference/best_reference.txt",
        ncbi_ref = f"{OUTPUT_DIR}/06_anomaly/reference/ncbi_reference.fasta",
        ref_gff = f"{OUTPUT_DIR}/06_anomaly/reference/ncbi_reference.gff"
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
        
        # Read reference path - prefer local, fallback to NCBI download
        with open(input.best_ref) as f:
            ref_path = f.read().strip()
        
        # If no local reference, use NCBI download
        if ref_path == "none" or not os.path.exists(ref_path):
            if os.path.exists(input.ncbi_ref) and os.path.getsize(input.ncbi_ref) > 0:
                ref_path = input.ncbi_ref
        
        alignment = {
            "sample_id": SAMPLE_ID,
            "reference": ref_path,
            "status": "no_reference",
            "aligned_regions": [],
            "unaligned_regions": [],
            "total_aligned_bp": 0,
            "insertions": [],
            "deletions": [],
            "deleted_genes": []
        }
        
        if ref_path != "none" and os.path.exists(ref_path) and os.path.getsize(ref_path) > 0:
            # Run nucmer
            cmd = f"nucmer --prefix={params.prefix} -t {threads} {ref_path} {input.assembly} 2>&1 | tee {log}"
            subprocess.run(cmd, shell=True)
            
            # Generate coords
            if os.path.exists(output.delta):
                cmd = f"show-coords -rcl {output.delta} > {output.coords} 2>> {log}"
                subprocess.run(cmd, shell=True)
                alignment["status"] = "aligned"
                
                # Parse coords to get aligned regions
                aligned_regions = []
                if os.path.exists(output.coords):
                    with open(output.coords) as f:
                        for line in f:
                            if line.startswith('=') or line.startswith('[') or '|' in line:
                                continue
                            parts = line.strip().split()
                            if len(parts) >= 12:
                                try:
                                    region = {
                                        "ref_start": int(parts[0]),
                                        "ref_end": int(parts[1]),
                                        "qry_start": int(parts[3]),
                                        "qry_end": int(parts[4]),
                                        "identity": float(parts[9]),
                                        "contig": parts[11]
                                    }
                                    aligned_regions.append(region)
                                    alignment["aligned_regions"].append(region)
                                    alignment["total_aligned_bp"] += abs(int(parts[1]) - int(parts[0]))
                                except (ValueError, IndexError):
                                    pass
                
                # Identify deletions (gaps in reference coverage)
                if aligned_regions:
                    # Sort by reference position
                    sorted_regions = sorted(aligned_regions, key=lambda x: x["ref_start"])
                    
                    # Find gaps between aligned regions
                    for i in range(len(sorted_regions) - 1):
                        gap_start = sorted_regions[i]["ref_end"]
                        gap_end = sorted_regions[i+1]["ref_start"]
                        gap_size = gap_end - gap_start
                        
                        if gap_size > 100:  # Minimum gap size
                            alignment["deletions"].append({
                                "ref_start": gap_start,
                                "ref_end": gap_end,
                                "size": gap_size
                            })
                
                # Map deletions to genes using reference GFF
                if os.path.exists(input.ref_gff) and os.path.getsize(input.ref_gff) > 0:
                    ref_genes = []
                    with open(input.ref_gff) as f:
                        for line in f:
                            if line.startswith('#'):
                                continue
                            parts = line.strip().split('\t')
                            if len(parts) >= 9 and parts[2] in ['gene', 'CDS']:
                                attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
                                ref_genes.append({
                                    "name": attrs.get("Name", attrs.get("locus_tag", attrs.get("ID", "unknown"))),
                                    "product": attrs.get("product", ""),
                                    "start": int(parts[3]),
                                    "end": int(parts[4])
                                })
                    
                    # Find genes in deleted regions
                    for deletion in alignment["deletions"]:
                        deleted_genes = []
                        for gene in ref_genes:
                            # Check if gene overlaps with deletion
                            if gene["start"] < deletion["ref_end"] and gene["end"] > deletion["ref_start"]:
                                deleted_genes.append({
                                    "name": gene["name"],
                                    "product": gene["product"]
                                })
                        
                        if deleted_genes:
                            alignment["deleted_genes"].append({
                                "deletion_start": deletion["ref_start"],
                                "deletion_end": deletion["ref_end"],
                                "genes": deleted_genes
                            })
        else:
            # Create empty files
            open(output.delta, 'w').close()
            open(output.coords, 'w').close()
        
        with open(output.alignment_json, 'w') as f:
            json.dump(alignment, f, indent=2)


# ============================================================================
# Synteny Analysis (Gene Order Comparison)
# ============================================================================

rule synteny_analysis:
    """Compare gene order between assembly and reference to detect rearrangements"""
    input:
        assembly_gff = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3",
        ref_gff = f"{OUTPUT_DIR}/06_anomaly/reference/ncbi_reference.gff",
        alignment = f"{OUTPUT_DIR}/06_anomaly/alignment/alignment.json"
    output:
        synteny_json = f"{OUTPUT_DIR}/06_anomaly/synteny/synteny_report.json"
    params:
        outdir = f"{OUTPUT_DIR}/06_anomaly/synteny"
    log:
        f"{OUTPUT_DIR}/00_logs/synteny.log"
    run:
        import json
        import os
        from collections import defaultdict
        
        os.makedirs(params.outdir, exist_ok=True)
        
        synteny = {
            "sample_id": SAMPLE_ID,
            "status": "not_analyzed",
            "total_genes_compared": 0,
            "genes_in_order": 0,
            "rearrangements": [],
            "inversions": [],
            "translocations": []
        }
        
        def parse_gff_genes(gff_file):
            """Extract gene names and positions from GFF"""
            genes = []
            if not os.path.exists(gff_file) or os.path.getsize(gff_file) == 0:
                return genes
            
            with open(gff_file) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 9 and parts[2] in ['gene', 'CDS']:
                        attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
                        gene_name = attrs.get("gene", attrs.get("Name", attrs.get("locus_tag", "")))
                        if gene_name:
                            genes.append({
                                "name": gene_name,
                                "contig": parts[0],
                                "start": int(parts[3]),
                                "end": int(parts[4]),
                                "strand": parts[6],
                                "product": attrs.get("product", "")
                            })
            return genes
        
        # Parse both GFF files
        assembly_genes = parse_gff_genes(input.assembly_gff)
        ref_genes = parse_gff_genes(input.ref_gff)
        
        if assembly_genes and ref_genes:
            synteny["status"] = "analyzed"
            
            # Create lookup by gene name
            ref_gene_order = {g["name"]: i for i, g in enumerate(ref_genes)}
            ref_gene_info = {g["name"]: g for g in ref_genes}
            
            # Compare gene order
            matched_genes = []
            for i, gene in enumerate(assembly_genes):
                if gene["name"] in ref_gene_order:
                    matched_genes.append({
                        "name": gene["name"],
                        "assembly_position": i,
                        "assembly_contig": gene["contig"],
                        "assembly_strand": gene["strand"],
                        "ref_position": ref_gene_order[gene["name"]],
                        "ref_strand": ref_gene_info[gene["name"]]["strand"]
                    })
            
            synteny["total_genes_compared"] = len(matched_genes)
            
            if len(matched_genes) >= 2:
                # Check for order violations
                prev_ref_pos = -1
                genes_in_order = 0
                
                for gene in matched_genes:
                    if gene["ref_position"] > prev_ref_pos:
                        genes_in_order += 1
                    else:
                        # Potential rearrangement
                        synteny["rearrangements"].append({
                            "gene": gene["name"],
                            "assembly_position": gene["assembly_position"],
                            "expected_ref_position": gene["ref_position"],
                            "type": "out_of_order"
                        })
                    
                    # Check for inversions (strand change)
                    if gene["assembly_strand"] != gene["ref_strand"]:
                        synteny["inversions"].append({
                            "gene": gene["name"],
                            "assembly_strand": gene["assembly_strand"],
                            "ref_strand": gene["ref_strand"]
                        })
                    
                    prev_ref_pos = gene["ref_position"]
                
                synteny["genes_in_order"] = genes_in_order
                
                # Detect translocations (genes from different ref regions on same contig)
                contig_genes = defaultdict(list)
                for gene in matched_genes:
                    contig_genes[gene["assembly_contig"]].append(gene)
                
                for contig, genes in contig_genes.items():
                    if len(genes) >= 2:
                        ref_positions = [g["ref_position"] for g in genes]
                        # Large jump in reference position = potential translocation
                        for i in range(len(ref_positions) - 1):
                            if abs(ref_positions[i+1] - ref_positions[i]) > 100:
                                synteny["translocations"].append({
                                    "contig": contig,
                                    "gene1": genes[i]["name"],
                                    "gene2": genes[i+1]["name"],
                                    "ref_distance": abs(ref_positions[i+1] - ref_positions[i])
                                })
        else:
            synteny["status"] = "missing_data"
            if not assembly_genes:
                synteny["error"] = "No assembly GFF available"
            elif not ref_genes:
                synteny["error"] = "No reference GFF available"
        
        with open(output.synteny_json, 'w') as f:
            json.dump(synteny, f, indent=2)


# ============================================================================
# Intergenic Region Analysis
# ============================================================================

rule extract_intergenic:
    """Extract intergenic regions for foreign element scanning"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        gff = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3"
    output:
        intergenic_fasta = f"{OUTPUT_DIR}/06_anomaly/intergenic/intergenic_regions.fasta"
    params:
        min_length = 100,  # Minimum intergenic region size
        max_length = 10000  # Maximum (avoid very long regions)
    run:
        import os
        from Bio import SeqIO
        
        os.makedirs(os.path.dirname(output.intergenic_fasta), exist_ok=True)
        
        # Load assembly sequences
        contigs = {rec.id: rec for rec in SeqIO.parse(input.assembly, "fasta")}
        
        # Parse GFF to get gene positions per contig
        contig_features = defaultdict(list)
        if os.path.exists(input.gff):
            with open(input.gff) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 5 and parts[2] in ['CDS', 'gene', 'tRNA', 'rRNA']:
                        contig_features[parts[0]].append({
                            "start": int(parts[3]),
                            "end": int(parts[4])
                        })
        
        # Extract intergenic regions
        intergenic_seqs = []
        for contig_id, features in contig_features.items():
            if contig_id not in contigs:
                continue
            
            contig_seq = contigs[contig_id]
            
            # Sort features by start position
            features = sorted(features, key=lambda x: x["start"])
            
            # Extract regions between features
            for i in range(len(features) - 1):
                start = features[i]["end"]
                end = features[i+1]["start"]
                length = end - start
                
                if params.min_length <= length <= params.max_length:
                    seq = str(contig_seq.seq[start:end])
                    intergenic_seqs.append({
                        "id": f"{contig_id}_intergenic_{start}_{end}",
                        "contig": contig_id,
                        "start": start,
                        "end": end,
                        "seq": seq
                    })
        
        # Write FASTA
        with open(output.intergenic_fasta, 'w') as f:
            for region in intergenic_seqs:
                f.write(f">{region['id']}\n{region['seq']}\n")


rule blast_intergenic:
    """BLAST intergenic regions to detect foreign regulatory elements"""
    input:
        intergenic = f"{OUTPUT_DIR}/06_anomaly/intergenic/intergenic_regions.fasta"
    output:
        blast_out = f"{OUTPUT_DIR}/06_anomaly/intergenic/blast_results.tsv",
        intergenic_json = f"{OUTPUT_DIR}/06_anomaly/intergenic/intergenic_taxonomy.json"
    params:
        blast_db = lambda w: config.get("databases", {}).get("nt", ""),
        evalue = 1e-5,
        max_targets = 3
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/blast_intergenic.log"
    conda:
        "../envs/anomaly.yaml"
    run:
        import json
        import subprocess
        import os
        
        intergenic_taxonomy = {
            "sample_id": SAMPLE_ID,
            "status": "not_run",
            "total_regions": 0,
            "regions_with_hits": 0,
            "foreign_elements": [],
            "host_genus": None
        }
        
        blast_db = params.blast_db
        
        if os.path.exists(input.intergenic) and os.path.getsize(input.intergenic) > 0:
            # Count regions
            with open(input.intergenic) as f:
                intergenic_taxonomy["total_regions"] = sum(1 for line in f if line.startswith('>'))
            
            if blast_db and os.path.exists(blast_db + ".nal"):
                # Run BLASTN
                cmd = f"""blastn \
                    -query {input.intergenic} \
                    -db {blast_db} \
                    -out {output.blast_out} \
                    -outfmt "6 qseqid sseqid pident length evalue staxids sscinames" \
                    -evalue {params.evalue} \
                    -max_target_seqs {params.max_targets} \
                    -num_threads {threads} 2>&1 | tee {log}"""
                subprocess.run(cmd, shell=True)
                
                intergenic_taxonomy["status"] = "completed"
                
                # Parse results and identify foreign elements
                region_hits = {}
                genus_counts = defaultdict(int)
                
                if os.path.exists(output.blast_out):
                    with open(output.blast_out) as f:
                        for line in f:
                            parts = line.strip().split('\t')
                            if len(parts) >= 7:
                                region_id = parts[0]
                                pident = float(parts[2])
                                organism = parts[6] if len(parts) > 6 else "Unknown"
                                genus = organism.split()[0] if organism else "Unknown"
                                
                                if region_id not in region_hits:
                                    region_hits[region_id] = {
                                        "organism": organism,
                                        "genus": genus,
                                        "identity": pident
                                    }
                                    genus_counts[genus] += 1
                
                intergenic_taxonomy["regions_with_hits"] = len(region_hits)
                
                # Determine host genus (from gene taxonomy if available)
                # For now, use most common genus
                if genus_counts:
                    host_genus = max(genus_counts, key=genus_counts.get)
                    intergenic_taxonomy["host_genus"] = host_genus
                    
                    # Flag foreign intergenic regions
                    for region_id, hit in region_hits.items():
                        if hit["genus"] != host_genus and hit["genus"] != "Unknown" and hit["identity"] > 80:
                            # Parse region coordinates from ID
                            parts = region_id.split('_')
                            intergenic_taxonomy["foreign_elements"].append({
                                "region_id": region_id,
                                "organism": hit["organism"],
                                "genus": hit["genus"],
                                "identity": hit["identity"],
                                "type": "foreign_regulatory_element"
                            })
            else:
                intergenic_taxonomy["status"] = "no_database"
                open(output.blast_out, 'w').close()
        else:
            intergenic_taxonomy["status"] = "no_intergenic_regions"
            open(output.blast_out, 'w').close()
        
        with open(output.intergenic_json, 'w') as f:
            json.dump(intergenic_taxonomy, f, indent=2)


# ============================================================================
# Integrated Mobile Element Detection
# ============================================================================

rule detect_integrated_elements:
    """Scan chromosome for plasmid backbone genes indicating integrated elements"""
    input:
        proteins = f"{OUTPUT_DIR}/06_anomaly/gene_taxonomy/proteins.faa",
        gff = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3"
    output:
        integrated_json = f"{OUTPUT_DIR}/06_anomaly/integrated_elements/integrated_elements.json"
    params:
        outdir = f"{OUTPUT_DIR}/06_anomaly/integrated_elements",
        # Plasmid backbone gene keywords
        plasmid_keywords = ["relaxase", "mobilization", "conjugal", "transfer", "replication", 
                          "repA", "repB", "repC", "traA", "traB", "traC", "traD", "traG",
                          "mobA", "mobB", "mobC", "oriT", "integrase", "transposase"]
    run:
        import json
        import os
        import re
        
        os.makedirs(params.outdir, exist_ok=True)
        
        integrated_elements = {
            "sample_id": SAMPLE_ID,
            "status": "analyzed",
            "plasmid_backbone_genes": [],
            "integrated_element_candidates": [],
            "total_backbone_genes": 0
        }
        
        # Parse GFF to get gene products and positions
        gene_info = {}
        if os.path.exists(input.gff):
            with open(input.gff) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 9 and parts[2] == 'CDS':
                        attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
                        gene_id = attrs.get("locus_tag", attrs.get("ID", ""))
                        product = attrs.get("product", "").lower()
                        
                        if gene_id:
                            gene_info[gene_id] = {
                                "contig": parts[0],
                                "start": int(parts[3]),
                                "end": int(parts[4]),
                                "strand": parts[6],
                                "product": product
                            }
        
        # Find plasmid backbone genes
        backbone_genes = []
        for gene_id, info in gene_info.items():
            for keyword in params.plasmid_keywords:
                if keyword.lower() in info["product"]:
                    backbone_genes.append({
                        "gene_id": gene_id,
                        "contig": info["contig"],
                        "start": info["start"],
                        "end": info["end"],
                        "product": info["product"],
                        "keyword_match": keyword
                    })
                    break
        
        integrated_elements["plasmid_backbone_genes"] = backbone_genes
        integrated_elements["total_backbone_genes"] = len(backbone_genes)
        
        # Cluster adjacent backbone genes to identify integrated elements
        if len(backbone_genes) >= 2:
            # Sort by contig and position
            sorted_genes = sorted(backbone_genes, key=lambda x: (x["contig"], x["start"]))
            
            # Find clusters (genes within 20kb of each other)
            clusters = []
            current_cluster = [sorted_genes[0]]
            
            for gene in sorted_genes[1:]:
                prev_gene = current_cluster[-1]
                
                if (gene["contig"] == prev_gene["contig"] and 
                    gene["start"] - prev_gene["end"] < 20000):
                    current_cluster.append(gene)
                else:
                    if len(current_cluster) >= 2:
                        clusters.append(current_cluster)
                    current_cluster = [gene]
            
            if len(current_cluster) >= 2:
                clusters.append(current_cluster)
            
            # Report clusters as integrated element candidates
            for i, cluster in enumerate(clusters):
                integrated_elements["integrated_element_candidates"].append({
                    "id": f"integrated_element_{i+1}",
                    "contig": cluster[0]["contig"],
                    "start": min(g["start"] for g in cluster),
                    "end": max(g["end"] for g in cluster),
                    "num_backbone_genes": len(cluster),
                    "genes": [g["gene_id"] for g in cluster],
                    "products": [g["product"] for g in cluster],
                    "type": "possible_integrated_plasmid_or_ICE"
                })
        
        with open(output.integrated_json, 'w') as f:
            json.dump(integrated_elements, f, indent=2)


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
        max_targets = 3,
        min_identity_threshold = 50  # Below this = novel/synthetic gene
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
            "novel_genes": [],  # Genes with no/weak BLAST hits
            "taxonomy_distribution": {},
            "total_genes": 0,
            "genes_with_hits": 0
        }
        
        blast_db = params.blast_db
        
        # Get all gene IDs from input proteins
        all_gene_ids = set()
        if os.path.exists(input.proteins) and os.path.getsize(input.proteins) > 0:
            with open(input.proteins) as f:
                for line in f:
                    if line.startswith('>'):
                        gene_id = line[1:].split()[0]
                        all_gene_ids.add(gene_id)
            gene_taxonomy["total_genes"] = len(all_gene_ids)
        
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
            genes_with_hits = set()
            
            with open(output.blast_out) as f:
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 7:
                        gene_id = parts[0]
                        pident = float(parts[2])
                        organism = parts[6] if len(parts) > 6 else "Unknown"
                        genus = organism.split()[0] if organism else "Unknown"
                        
                        if gene_id not in gene_taxonomy["genes"]:
                            genes_with_hits.add(gene_id)
                            gene_taxonomy["genes"][gene_id] = {
                                "organism": organism,
                                "genus": genus,
                                "identity": pident,
                                "status": "matched"
                            }
                            genus_counts[genus] += 1
            
            gene_taxonomy["genes_with_hits"] = len(genes_with_hits)
            
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
                    
                    # Flag low-identity hits as potentially novel
                    if data["identity"] < params.min_identity_threshold:
                        gene_taxonomy["novel_genes"].append({
                            "gene_id": gene_id,
                            "best_hit": data["organism"],
                            "identity": data["identity"],
                            "reason": "low_identity_hit"
                        })
            
            # Find genes with NO BLAST hits at all - these are suspicious
            genes_without_hits = all_gene_ids - genes_with_hits
            for gene_id in genes_without_hits:
                gene_taxonomy["genes"][gene_id] = {
                    "organism": None,
                    "genus": None,
                    "identity": 0,
                    "status": "no_hit"
                }
                gene_taxonomy["novel_genes"].append({
                    "gene_id": gene_id,
                    "best_hit": None,
                    "identity": 0,
                    "reason": "no_blast_hit"
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
    """Combine all anomaly detection results and cluster foreign genes"""
    input:
        ref_info = f"{OUTPUT_DIR}/06_anomaly/reference/reference_info.json",
        alignment = f"{OUTPUT_DIR}/06_anomaly/alignment/alignment.json",
        gene_taxonomy = f"{OUTPUT_DIR}/06_anomaly/gene_taxonomy/gene_taxonomy.json",
        synteny = f"{OUTPUT_DIR}/06_anomaly/synteny/synteny_report.json",
        intergenic = f"{OUTPUT_DIR}/06_anomaly/intergenic/intergenic_taxonomy.json",
        integrated = f"{OUTPUT_DIR}/06_anomaly/integrated_elements/integrated_elements.json",
        gff = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3"
    output:
        anomaly_summary = f"{OUTPUT_DIR}/06_anomaly/anomaly_summary.json"
    params:
        min_cluster_size = 3  # Minimum adjacent foreign genes to flag as integrated element
    run:
        import json
        import os
        
        # Helper function to load JSON safely
        def load_json_safe(path, default={}):
            if os.path.exists(path):
                try:
                    with open(path) as f:
                        return json.load(f)
                except:
                    return default
            return default
        
        # Load all inputs
        ref_info = load_json_safe(input.ref_info)
        alignment = load_json_safe(input.alignment)
        gene_taxonomy = load_json_safe(input.gene_taxonomy)
        synteny = load_json_safe(input.synteny)
        intergenic = load_json_safe(input.intergenic)
        integrated = load_json_safe(input.integrated)
        
        # Parse GFF to get gene positions for clustering
        gene_positions = {}
        if os.path.exists(input.gff):
            with open(input.gff) as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 9 and parts[2] == 'CDS':
                        attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
                        gene_id = attrs.get("locus_tag", attrs.get("ID", ""))
                        if gene_id:
                            gene_positions[gene_id] = {
                                "contig": parts[0],
                                "start": int(parts[3]),
                                "end": int(parts[4])
                            }
        
        # Cluster foreign genes by position
        foreign_gene_clusters = []
        foreign_genes = gene_taxonomy.get("foreign_genes", [])
        
        if foreign_genes and gene_positions:
            # Add position info to foreign genes
            positioned_foreign = []
            for fg in foreign_genes:
                gene_id = fg["gene_id"]
                if gene_id in gene_positions:
                    positioned_foreign.append({
                        **fg,
                        "contig": gene_positions[gene_id]["contig"],
                        "start": gene_positions[gene_id]["start"],
                        "end": gene_positions[gene_id]["end"]
                    })
            
            # Sort by contig and position
            positioned_foreign.sort(key=lambda x: (x["contig"], x["start"]))
            
            # Cluster adjacent foreign genes (within 10kb)
            if positioned_foreign:
                current_cluster = [positioned_foreign[0]]
                
                for gene in positioned_foreign[1:]:
                    prev = current_cluster[-1]
                    
                    if (gene["contig"] == prev["contig"] and 
                        gene["start"] - prev["end"] < 10000):
                        current_cluster.append(gene)
                    else:
                        if len(current_cluster) >= params.min_cluster_size:
                            # Check if genes are from same genus
                            genera = [g["genus"] for g in current_cluster]
                            dominant_genus = max(set(genera), key=genera.count)
                            
                            foreign_gene_clusters.append({
                                "contig": current_cluster[0]["contig"],
                                "start": min(g["start"] for g in current_cluster),
                                "end": max(g["end"] for g in current_cluster),
                                "num_genes": len(current_cluster),
                                "genes": [g["gene_id"] for g in current_cluster],
                                "dominant_genus": dominant_genus,
                                "type": "integrated_foreign_element"
                            })
                        current_cluster = [gene]
                
                # Don't forget last cluster
                if len(current_cluster) >= params.min_cluster_size:
                    genera = [g["genus"] for g in current_cluster]
                    dominant_genus = max(set(genera), key=genera.count)
                    foreign_gene_clusters.append({
                        "contig": current_cluster[0]["contig"],
                        "start": min(g["start"] for g in current_cluster),
                        "end": max(g["end"] for g in current_cluster),
                        "num_genes": len(current_cluster),
                        "genes": [g["gene_id"] for g in current_cluster],
                        "dominant_genus": dominant_genus,
                        "type": "integrated_foreign_element"
                    })
        
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
                "deletions": alignment.get("deletions", []),
                "deleted_genes": alignment.get("deleted_genes", [])
            },
            
            # Gene taxonomy
            "gene_taxonomy": {
                "status": gene_taxonomy.get("status"),
                "host_genus": gene_taxonomy.get("host_genus"),
                "foreign_genes": gene_taxonomy.get("foreign_genes", []),
                "num_foreign_genes": len(gene_taxonomy.get("foreign_genes", [])),
                "novel_genes": gene_taxonomy.get("novel_genes", []),
                "num_novel_genes": len(gene_taxonomy.get("novel_genes", [])),
                "taxonomy_distribution": gene_taxonomy.get("taxonomy_distribution", {})
            },
            
            # Foreign gene clusters
            "foreign_gene_clusters": foreign_gene_clusters,
            
            # Synteny analysis
            "synteny": {
                "status": synteny.get("status"),
                "total_genes_compared": synteny.get("total_genes_compared", 0),
                "genes_in_order": synteny.get("genes_in_order", 0),
                "rearrangements": synteny.get("rearrangements", []),
                "inversions": synteny.get("inversions", []),
                "translocations": synteny.get("translocations", [])
            },
            
            # Intergenic region analysis
            "intergenic": {
                "status": intergenic.get("status"),
                "foreign_elements": intergenic.get("foreign_elements", []),
                "num_foreign_elements": len(intergenic.get("foreign_elements", []))
            },
            
            # Integrated mobile elements
            "integrated_elements": {
                "candidates": integrated.get("integrated_element_candidates", []),
                "num_candidates": len(integrated.get("integrated_element_candidates", []))
            },
            
            # Anomaly flags
            "anomaly_detected": False,
            "anomaly_types": [],
            "severity": "none"
        }
        
        # Check for anomalies and determine severity
        anomaly_score = 0
        
        if summary["gene_taxonomy"]["num_foreign_genes"] > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("foreign_gene_insertion")
            anomaly_score += summary["gene_taxonomy"]["num_foreign_genes"]
        
        if summary["gene_taxonomy"]["num_novel_genes"] > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("novel_synthetic_genes")
            anomaly_score += summary["gene_taxonomy"]["num_novel_genes"] * 2  # Weight novel genes higher
        
        if len(foreign_gene_clusters) > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("integrated_foreign_element")
            anomaly_score += len(foreign_gene_clusters) * 3  # Weight clusters high
        
        if len(summary["alignment"]["deletions"]) > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("structural_deletion")
            anomaly_score += len(summary["alignment"]["deletions"])
        
        if len(summary["alignment"]["deleted_genes"]) > 0:
            summary["anomaly_types"].append("gene_deletion")
        
        if len(synteny.get("rearrangements", [])) > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("gene_rearrangement")
            anomaly_score += 1
        
        if len(synteny.get("inversions", [])) > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("gene_inversion")
            anomaly_score += 1
        
        if summary["intergenic"]["num_foreign_elements"] > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("foreign_regulatory_element")
            anomaly_score += summary["intergenic"]["num_foreign_elements"]
        
        if summary["integrated_elements"]["num_candidates"] > 0:
            summary["anomaly_detected"] = True
            summary["anomaly_types"].append("possible_integrated_plasmid")
            anomaly_score += summary["integrated_elements"]["num_candidates"] * 2
        
        # Determine severity
        if anomaly_score == 0:
            summary["severity"] = "none"
        elif anomaly_score <= 3:
            summary["severity"] = "low"
        elif anomaly_score <= 10:
            summary["severity"] = "medium"
        else:
            summary["severity"] = "high"
        
        with open(output.anomaly_summary, 'w') as f:
            json.dump(summary, f, indent=2)
