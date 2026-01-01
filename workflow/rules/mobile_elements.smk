"""
============================================================================
Module 5d: Mobile Genetic Elements (MGE)
============================================================================
Detect plasmids and mobile elements using MOB-suite
"""

# ============================================================================
# MOB-suite
# ============================================================================

rule mob_recon:
    """Identify and type plasmids with MOB-suite"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
    output:
        contig_report = f"{OUTPUT_DIR}/05_specialized/mge/contig_report.txt",
        mobtyper = f"{OUTPUT_DIR}/05_specialized/mge/mobtyper_results.txt",
        chromosome = f"{OUTPUT_DIR}/05_specialized/mge/chromosome.fasta",
        flag = f"{OUTPUT_DIR}/05_specialized/mge/mob_suite_complete.flag",
    params:
        outdir = f"{OUTPUT_DIR}/05_specialized/mge",
        min_length = lambda w: config.get("mobile_elements", {}).get("mob_suite", {}).get("min_length", 1000),
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/mob_recon.log"
    conda:
        "../envs/specialized.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run MOB-recon
        mob_recon \
            --infile {input.assembly} \
            --outdir {params.outdir} \
            --num_threads {threads} \
            --force \
            2>&1 | tee {log}
        
        # Ensure output files exist
        touch {output.contig_report} {output.mobtyper} {output.chromosome}
        touch {output.flag}
        """


rule parse_mob_results:
    """Parse MOB-suite results into structured format"""
    input:
        contig_report = rules.mob_recon.output.contig_report,
        mobtyper = rules.mob_recon.output.mobtyper,
        flag = rules.mob_recon.output.flag,
    output:
        json = f"{OUTPUT_DIR}/05_specialized/mge/mge_summary.json",
    run:
        import json
        import csv
        import glob
        
        mge_data = {
            "sample_id": SAMPLE_ID,
            "tool": "MOB-suite",
            "plasmids": [],
            "summary": {
                "total_plasmids": 0,
                "chromosome_contigs": 0,
                "plasmid_contigs": 0,
            }
        }
        
        # Parse contig report
        if os.path.exists(input.contig_report) and os.path.getsize(input.contig_report) > 0:
            with open(input.contig_report) as f:
                reader = csv.DictReader(f, delimiter='\t')
                chromosome_contigs = 0
                plasmid_contigs = 0
                
                for row in reader:
                    if row.get("cluster_id", "").lower() == "chromosome":
                        chromosome_contigs += 1
                    else:
                        plasmid_contigs += 1
                
                mge_data["summary"]["chromosome_contigs"] = chromosome_contigs
                mge_data["summary"]["plasmid_contigs"] = plasmid_contigs
        
        # Parse mobtyper results
        if os.path.exists(input.mobtyper) and os.path.getsize(input.mobtyper) > 0:
            with open(input.mobtyper) as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    plasmid = {
                        "cluster_id": row.get("cluster_id", ""),
                        "size": int(row.get("size", 0) or 0),
                        "gc_content": float(row.get("gc", 0) or 0),
                        "rep_type": row.get("rep_type(s)", ""),
                        "relaxase_type": row.get("relaxase_type(s)", ""),
                        "mob_type": row.get("mpf_type", ""),
                        "predicted_mobility": row.get("predicted_mobility", ""),
                        "primary_cluster_id": row.get("primary_cluster_id", ""),
                    }
                    mge_data["plasmids"].append(plasmid)
        
        # Check for plasmid FASTA files
        plasmid_fastas = glob.glob(f"{OUTPUT_DIR}/05_specialized/mge/plasmid_*.fasta")
        mge_data["summary"]["total_plasmids"] = len(plasmid_fastas)
        
        with open(output.json, 'w') as f:
            json.dump(mge_data, f, indent=2)

