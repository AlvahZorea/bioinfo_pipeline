"""
============================================================================
Module 5b: Virulence Factors
============================================================================
Detect virulence factors using BLAST against VFDB
"""

# ============================================================================
# VFDB BLAST
# ============================================================================

rule vfdb_blast:
    """BLAST proteins against VFDB for virulence factors"""
    input:
        proteins = f"{OUTPUT_DIR}/04_annotation/proteins.faa",
    output:
        results = f"{OUTPUT_DIR}/05_specialized/virulence/vfdb_results.tsv",
    params:
        db = lambda w: get_database_path(config, "vfdb", "protein"),
        identity = lambda w: config.get("virulence", {}).get("vfdb", {}).get("identity", 80),
        coverage = lambda w: config.get("virulence", {}).get("vfdb", {}).get("coverage", 60),
        evalue = lambda w: config.get("virulence", {}).get("vfdb", {}).get("evalue", 1e-10),
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/vfdb_blast.log"
    conda:
        "../envs/specialized.yaml"
    shell:
        """
        mkdir -p $(dirname {output.results})
        
        # Check if input has sequences
        if [ -s {input.proteins} ]; then
            blastp \
                -query {input.proteins} \
                -db {params.db} \
                -out {output.results} \
                -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs stitle" \
                -evalue {params.evalue} \
                -num_threads {threads} \
                -max_target_seqs 5 \
                2>&1 | tee {log}
        else
            echo "No proteins to search" | tee {log}
            touch {output.results}
        fi
        """


rule filter_vfdb_results:
    """Filter VFDB BLAST results by identity and coverage"""
    input:
        rules.vfdb_blast.output.results
    output:
        filtered = f"{OUTPUT_DIR}/05_specialized/virulence/vfdb_filtered.tsv",
    params:
        identity = lambda w: config.get("virulence", {}).get("vfdb", {}).get("identity", 80),
        coverage = lambda w: config.get("virulence", {}).get("vfdb", {}).get("coverage", 60),
    run:
        with open(input[0]) as f_in, open(output.filtered, 'w') as f_out:
            f_out.write("query_id\tsubject_id\tidentity\talign_length\tmismatches\tgaps\tq_start\tq_end\ts_start\ts_end\tevalue\tbitscore\tcoverage\tdescription\n")
            for line in f_in:
                parts = line.strip().split('\t')
                if len(parts) >= 13:
                    identity = float(parts[2])
                    coverage = float(parts[12])
                    if identity >= params.identity and coverage >= params.coverage:
                        f_out.write(line)


rule parse_virulence_results:
    """Parse virulence results into structured format"""
    input:
        filtered = rules.filter_vfdb_results.output.filtered,
    output:
        json = f"{OUTPUT_DIR}/05_specialized/virulence/virulence_summary.json",
    run:
        import json
        import re
        
        virulence_data = {
            "sample_id": SAMPLE_ID,
            "tool": "BLAST vs VFDB",
            "virulence_factors": [],
            "categories": {},
            "summary": {
                "total_factors": 0,
            }
        }
        
        seen = set()  # Track unique factors
        
        if os.path.exists(input.filtered):
            with open(input.filtered) as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 14:
                        query_id = parts[0]
                        subject_id = parts[1]
                        identity = float(parts[2])
                        coverage = float(parts[12])
                        description = parts[13] if len(parts) > 13 else ""
                        
                        # Extract VF category from description if possible
                        category = "Unknown"
                        if "[" in description and "]" in description:
                            match = re.search(r'\[([^\]]+)\]', description)
                            if match:
                                category = match.group(1)
                        
                        # Only add unique factors (by subject_id)
                        if subject_id not in seen:
                            seen.add(subject_id)
                            
                            factor = {
                                "query": query_id,
                                "vfdb_id": subject_id,
                                "identity": identity,
                                "coverage": coverage,
                                "description": description,
                                "category": category,
                            }
                            virulence_data["virulence_factors"].append(factor)
                            
                            # Count by category
                            virulence_data["categories"][category] = virulence_data["categories"].get(category, 0) + 1
        
        virulence_data["summary"]["total_factors"] = len(virulence_data["virulence_factors"])
        
        with open(output.json, 'w') as f:
            json.dump(virulence_data, f, indent=2)

