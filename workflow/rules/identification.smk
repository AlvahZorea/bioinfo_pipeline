"""
============================================================================
Module 3: Identification
============================================================================
Taxonomic identification and closest reference detection

Steps:
    1. Kraken2 classification on assembly
    2. fastANI for closest reference
    3. geNomad for virus/plasmid identification
    4. CheckV for viral genome quality
    5. Determine organism type (bacteria vs virus)
"""

# ============================================================================
# Kraken2 Classification
# ============================================================================

rule kraken2_assembly:
    """Classify assembled contigs with Kraken2"""
    input:
        f"{OUTPUT_DIR}/02_assembly/contigs.fasta"
    output:
        report = f"{OUTPUT_DIR}/03_identification/kraken2_report.txt",
        output = f"{OUTPUT_DIR}/03_identification/kraken2_output.txt",
    params:
        db = lambda w: get_database_path(config, "kraken2", "standard"),
        confidence = lambda w: config.get("identification", {}).get("kraken2", {}).get("confidence", 0.1),
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/kraken2_assembly.log"
    conda:
        "../envs/identification.yaml"
    shell:
        """
        mkdir -p $(dirname {output.report})
        
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --confidence {params.confidence} \
            --report {output.report} \
            --output {output.output} \
            {input} \
            2>&1 | tee {log}
        """


rule parse_kraken_taxonomy:
    """Parse Kraken2 results to extract taxonomy with minority species detection"""
    input:
        report = rules.kraken2_assembly.output.report,
        output = rules.kraken2_assembly.output.output,
    output:
        taxonomy = f"{OUTPUT_DIR}/03_identification/kraken2_taxonomy.json",
    params:
        minority_threshold = 0.5  # Percentage threshold for flagging minority species
    run:
        import json
        
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
            "minority_species_alert": False,
            "minority_species": []
        }
        
        rank_map = {'D': 'domain', 'P': 'phylum', 'C': 'class', 'O': 'order', 
                    'F': 'family', 'G': 'genus', 'S': 'species'}
        
        # Collect all genus-level hits for minority detection
        all_genera = []
        
        with open(input.report) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    pct = float(parts[0])
                    rank = parts[3]
                    taxid = parts[4]
                    name = parts[5].strip()
                    
                    # Store top hits
                    if pct > 1.0 and rank in ['G', 'S']:
                        taxonomy["top_hits"].append({
                            "name": name,
                            "taxid": taxid,
                            "rank": rank,
                            "percentage": pct
                        })
                    
                    # Collect all genus-level hits for minority detection
                    if rank == 'G' and pct >= params.minority_threshold:
                        all_genera.append({
                            "name": name,
                            "taxid": taxid,
                            "percentage": pct
                        })
                    
                    # Store highest percentage for each rank
                    if rank in rank_map and pct > 0:
                        rank_name = rank_map[rank]
                        if taxonomy[rank_name] is None or pct > taxonomy.get(f"{rank_name}_pct", 0):
                            taxonomy[rank_name] = name
                            taxonomy[f"{rank_name}_pct"] = pct
        
        # Sort top hits by percentage
        taxonomy["top_hits"] = sorted(taxonomy["top_hits"], key=lambda x: x["percentage"], reverse=True)[:10]
        
        # Detect minority species (any genus above threshold that isn't the dominant one)
        if len(all_genera) >= 2:
            # Sort by percentage
            all_genera = sorted(all_genera, key=lambda x: x["percentage"], reverse=True)
            dominant_genus = all_genera[0]["name"]
            
            # Flag non-dominant genera above threshold
            for genus in all_genera[1:]:
                if genus["percentage"] >= params.minority_threshold:
                    taxonomy["minority_species_alert"] = True
                    taxonomy["minority_species"].append({
                        "name": genus["name"],
                        "taxid": genus["taxid"],
                        "percentage": genus["percentage"],
                        "warning": f"Non-dominant genus at {genus['percentage']:.1f}% - may indicate contamination or chimeric genome"
                    })
        
        with open(output.taxonomy, 'w') as f:
            json.dump(taxonomy, f, indent=2)


# ============================================================================
# fastANI
# ============================================================================

rule fastani:
    """Calculate ANI to reference genomes"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
    output:
        result = f"{OUTPUT_DIR}/03_identification/fastani_results.txt",
    params:
        reference_dir = lambda w: get_database_path(config, "references"),
        min_fraction = lambda w: config.get("identification", {}).get("fastani", {}).get("min_fraction", 0.2),
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/fastani.log"
    conda:
        "../envs/identification.yaml"
    shell:
        """
        mkdir -p $(dirname {output.result})
        
        # Create reference list if directory provided
        if [ -d "{params.reference_dir}" ]; then
            find {params.reference_dir} -name "*.fasta" -o -name "*.fna" -o -name "*.fa" > /tmp/ref_list.txt
            
            if [ -s /tmp/ref_list.txt ]; then
                fastANI \
                    -q {input.assembly} \
                    --rl /tmp/ref_list.txt \
                    -o {output.result} \
                    -t {threads} \
                    --minFraction {params.min_fraction} \
                    2>&1 | tee {log}
            else
                echo "No reference genomes found" > {output.result}
            fi
        else
            echo "Reference directory not configured" > {output.result}
        fi
        """


rule parse_fastani:
    """Parse fastANI results"""
    input:
        rules.fastani.output.result
    output:
        f"{OUTPUT_DIR}/03_identification/fastani_parsed.json"
    run:
        import json
        
        results = {
            "sample_id": SAMPLE_ID,
            "tool": "fastANI",
            "closest_references": []
        }
        
        with open(input[0]) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 5:
                    results["closest_references"].append({
                        "query": parts[0],
                        "reference": parts[1],
                        "ani": float(parts[2]),
                        "matched_fragments": int(parts[3]),
                        "total_fragments": int(parts[4])
                    })
        
        # Sort by ANI
        results["closest_references"] = sorted(
            results["closest_references"], 
            key=lambda x: x["ani"], 
            reverse=True
        )[:10]
        
        with open(output[0], 'w') as f:
            json.dump(results, f, indent=2)


# ============================================================================
# geNomad (Virus/Plasmid Detection)
# ============================================================================

rule genomad:
    """Identify viral and plasmid sequences with geNomad"""
    input:
        f"{OUTPUT_DIR}/02_assembly/contigs.fasta"
    output:
        summary = f"{OUTPUT_DIR}/03_identification/genomad/{SAMPLE_ID}_summary.json",
        virus_fasta = f"{OUTPUT_DIR}/03_identification/genomad/{SAMPLE_ID}_virus.fna",
        plasmid_fasta = f"{OUTPUT_DIR}/03_identification/genomad/{SAMPLE_ID}_plasmid.fna",
        flag = f"{OUTPUT_DIR}/03_identification/genomad/genomad_complete.flag",
    params:
        db = lambda w: get_database_path(config, "genomad"),
        outdir = f"{OUTPUT_DIR}/03_identification/genomad",
        min_score = lambda w: config.get("identification", {}).get("genomad", {}).get("min_score", 0.7),
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/genomad.log"
    conda:
        "../envs/identification.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        genomad end-to-end \
            {input} \
            {params.outdir} \
            {params.db} \
            --threads {threads} \
            --min-score {params.min_score} \
            2>&1 | tee {log}
        
        # Collect outputs (handle if files don't exist)
        touch {output.virus_fasta} {output.plasmid_fasta}
        
        # Find and copy virus sequences
        if ls {params.outdir}/*_virus/*_virus.fna 1> /dev/null 2>&1; then
            cat {params.outdir}/*_virus/*_virus.fna > {output.virus_fasta}
        fi
        
        # Find and copy plasmid sequences
        if ls {params.outdir}/*_plasmid/*_plasmid.fna 1> /dev/null 2>&1; then
            cat {params.outdir}/*_plasmid/*_plasmid.fna > {output.plasmid_fasta}
        fi
        
        # Create summary JSON
        echo '{{"sample_id": "{SAMPLE_ID}", "tool": "geNomad"}}' > {output.summary}
        
        touch {output.flag}
        """


# ============================================================================
# CheckV (Viral Quality)
# ============================================================================

rule checkv:
    """Assess viral genome quality with CheckV"""
    input:
        virus_fasta = rules.genomad.output.virus_fasta,
    output:
        quality = f"{OUTPUT_DIR}/03_identification/checkv/quality_summary.tsv",
        completeness = f"{OUTPUT_DIR}/03_identification/checkv/completeness.tsv",
        flag = f"{OUTPUT_DIR}/03_identification/checkv/checkv_complete.flag",
    params:
        db = lambda w: get_database_path(config, "checkv"),
        outdir = f"{OUTPUT_DIR}/03_identification/checkv",
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/checkv.log"
    conda:
        "../envs/identification.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Only run if there are viral sequences
        if [ -s {input.virus_fasta} ]; then
            checkv end_to_end \
                {input.virus_fasta} \
                {params.outdir} \
                -d {params.db} \
                -t {threads} \
                2>&1 | tee {log}
        else
            echo "No viral sequences to analyze" | tee {log}
            echo -e "contig_id\\tquality" > {output.quality}
            echo -e "contig_id\\tcompleteness" > {output.completeness}
        fi
        
        touch {output.quality} {output.completeness} {output.flag}
        """


# ============================================================================
# Organism Type Determination
# ============================================================================

rule determine_organism_type:
    """Determine if sample is bacteria, virus, or mixed"""
    input:
        kraken_taxonomy = rules.parse_kraken_taxonomy.output.taxonomy,
        genomad_flag = rules.genomad.output.flag,
        checkv_flag = rules.checkv.output.flag,
    output:
        organism_type = f"{OUTPUT_DIR}/03_identification/organism_type.txt",
        taxonomy = f"{OUTPUT_DIR}/03_identification/taxonomy.json",
    run:
        import json
        import os
        
        # Load Kraken taxonomy
        with open(input.kraken_taxonomy) as f:
            kraken_tax = json.load(f)
        
        # Check for viral content via geNomad
        virus_fasta = f"{OUTPUT_DIR}/03_identification/genomad/{SAMPLE_ID}_virus.fna"
        has_virus = os.path.exists(virus_fasta) and os.path.getsize(virus_fasta) > 0
        
        # Determine organism type
        domain = kraken_tax.get("domain", "").lower()
        
        if "virus" in domain or has_virus:
            organism_type = "virus"
        elif "bacteria" in domain:
            organism_type = "bacteria"
        elif "archaea" in domain:
            organism_type = "archaea"
        else:
            # Default to bacteria if unclear
            organism_type = "bacteria"
        
        # Write organism type
        with open(output.organism_type, 'w') as f:
            f.write(organism_type)
        
        # Create combined taxonomy JSON
        taxonomy = {
            "sample_id": SAMPLE_ID,
            "organism_type": organism_type,
            "kraken2": kraken_tax,
            "has_viral_content": has_virus,
        }
        
        with open(output.taxonomy, 'w') as f:
            json.dump(taxonomy, f, indent=2)

