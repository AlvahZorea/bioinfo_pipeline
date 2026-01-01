"""
============================================================================
Module 5a: Antimicrobial Resistance (AMR)
============================================================================
Detect AMR genes and point mutations using AMRFinderPlus
"""

# ============================================================================
# AMRFinderPlus
# ============================================================================

rule amrfinderplus:
    """Detect AMR genes with AMRFinderPlus"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        proteins = f"{OUTPUT_DIR}/04_annotation/proteins.faa",
        gff = f"{OUTPUT_DIR}/04_annotation/annotation.gff",
    output:
        results = f"{OUTPUT_DIR}/05_specialized/amr/amrfinder_results.tsv",
        mutations = f"{OUTPUT_DIR}/05_specialized/amr/amrfinder_mutations.tsv",
    params:
        db = lambda w: get_database_path(config, "amrfinder"),
        ident_min = lambda w: config.get("amr", {}).get("amrfinderplus", {}).get("ident_min", 0.9),
        coverage_min = lambda w: config.get("amr", {}).get("amrfinderplus", {}).get("coverage_min", 0.5),
        organism = lambda w: config.get("amr", {}).get("amrfinderplus", {}).get("organism"),
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/amrfinderplus.log"
    conda:
        "../envs/specialized.yaml"
    shell:
        """
        mkdir -p $(dirname {output.results})
        
        # Build organism flag if specified
        ORGANISM_FLAG=""
        if [ -n "{params.organism}" ] && [ "{params.organism}" != "None" ]; then
            ORGANISM_FLAG="--organism {params.organism}"
        fi
        
        # Run AMRFinderPlus
        amrfinder \
            --nucleotide {input.assembly} \
            --protein {input.proteins} \
            --gff {input.gff} \
            --database {params.db} \
            --threads {threads} \
            --ident_min {params.ident_min} \
            --coverage_min {params.coverage_min} \
            --output {output.results} \
            --mutation_all {output.mutations} \
            --plus \
            $ORGANISM_FLAG \
            2>&1 | tee {log}
        
        # Ensure output files exist
        touch {output.results} {output.mutations}
        """


rule parse_amr_results:
    """Parse AMRFinderPlus results into structured format"""
    input:
        results = rules.amrfinderplus.output.results,
        mutations = rules.amrfinderplus.output.mutations,
    output:
        json = f"{OUTPUT_DIR}/05_specialized/amr/amr_summary.json",
    run:
        import json
        import csv
        
        amr_data = {
            "sample_id": SAMPLE_ID,
            "tool": "AMRFinderPlus",
            "genes": [],
            "mutations": [],
            "drug_classes": set(),
            "summary": {
                "total_genes": 0,
                "total_mutations": 0,
            }
        }
        
        # Parse gene results
        if os.path.exists(input.results) and os.path.getsize(input.results) > 0:
            with open(input.results) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    gene_info = {
                        "gene_symbol": row.get("Gene symbol", ""),
                        "sequence_name": row.get("Sequence name", ""),
                        "element_type": row.get("Element type", ""),
                        "element_subtype": row.get("Element subtype", ""),
                        "drug_class": row.get("Class", ""),
                        "subclass": row.get("Subclass", ""),
                        "identity": float(row.get("% Identity of coverage", 0) or 0),
                        "coverage": float(row.get("% Coverage of reference", 0) or 0),
                        "contig": row.get("Contig id", ""),
                        "start": row.get("Start", ""),
                        "stop": row.get("Stop", ""),
                    }
                    amr_data["genes"].append(gene_info)
                    if gene_info["drug_class"]:
                        amr_data["drug_classes"].add(gene_info["drug_class"])
        
        # Parse mutation results
        if os.path.exists(input.mutations) and os.path.getsize(input.mutations) > 0:
            with open(input.mutations) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    mutation_info = {
                        "gene": row.get("Gene symbol", ""),
                        "mutation": row.get("Point mutation", ""),
                        "drug_class": row.get("Class", ""),
                    }
                    if mutation_info["mutation"]:
                        amr_data["mutations"].append(mutation_info)
        
        # Update summary
        amr_data["summary"]["total_genes"] = len(amr_data["genes"])
        amr_data["summary"]["total_mutations"] = len(amr_data["mutations"])
        amr_data["drug_classes"] = list(amr_data["drug_classes"])
        
        with open(output.json, 'w') as f:
            json.dump(amr_data, f, indent=2)

