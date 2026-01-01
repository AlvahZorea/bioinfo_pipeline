"""
============================================================================
Module 5e: Prophage Detection
============================================================================
Detect prophage regions using PhiSpy
"""

# ============================================================================
# PhiSpy
# ============================================================================

rule phispy:
    """Detect prophages with PhiSpy"""
    input:
        gbk = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gbff",
    output:
        prophages = f"{OUTPUT_DIR}/05_specialized/prophages/prophage_coordinates.tsv",
        fasta = f"{OUTPUT_DIR}/05_specialized/prophages/prophage.fasta",
        flag = f"{OUTPUT_DIR}/05_specialized/prophages/phispy_complete.flag",
    params:
        outdir = f"{OUTPUT_DIR}/05_specialized/prophages",
        min_contig = lambda w: config.get("prophages", {}).get("phispy", {}).get("min_contig_length", 5000),
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/phispy.log"
    conda:
        "../envs/specialized.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Run PhiSpy
        PhiSpy.py \
            {input.gbk} \
            -o {params.outdir} \
            --threads {threads} \
            --output_choice 4 \
            2>&1 | tee {log} || true
        
        # Handle outputs (PhiSpy may not find any prophages)
        if [ -f "{params.outdir}/prophage_coordinates.tsv" ]; then
            cp {params.outdir}/prophage_coordinates.tsv {output.prophages}
        else
            echo -e "pp_no\\tcontig\\tstart\\tend\\tattL\\tattR\\tgc\\tstatus\\tsequence" > {output.prophages}
        fi
        
        if [ -f "{params.outdir}/prophage.fasta" ]; then
            cp {params.outdir}/prophage.fasta {output.fasta}
        else
            touch {output.fasta}
        fi
        
        touch {output.flag}
        """


rule parse_prophage_results:
    """Parse PhiSpy results into structured format"""
    input:
        prophages = rules.phispy.output.prophages,
        fasta = rules.phispy.output.fasta,
        flag = rules.phispy.output.flag,
    output:
        json = f"{OUTPUT_DIR}/05_specialized/prophages/prophage_summary.json",
    run:
        import json
        from Bio import SeqIO
        
        prophage_data = {
            "sample_id": SAMPLE_ID,
            "tool": "PhiSpy",
            "prophages": [],
            "summary": {
                "total_prophages": 0,
                "total_length": 0,
            }
        }
        
        # Parse prophage coordinates
        if os.path.exists(input.prophages) and os.path.getsize(input.prophages) > 0:
            with open(input.prophages) as f:
                header = f.readline()  # Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        try:
                            prophage = {
                                "id": parts[0],
                                "contig": parts[1],
                                "start": int(parts[2]),
                                "end": int(parts[3]),
                                "length": int(parts[3]) - int(parts[2]),
                                "attL": parts[4] if len(parts) > 4 else "",
                                "attR": parts[5] if len(parts) > 5 else "",
                                "gc_content": float(parts[6]) if len(parts) > 6 and parts[6] else 0,
                                "status": parts[7] if len(parts) > 7 else "",
                            }
                            prophage_data["prophages"].append(prophage)
                            prophage_data["summary"]["total_length"] += prophage["length"]
                        except (ValueError, IndexError):
                            continue
        
        prophage_data["summary"]["total_prophages"] = len(prophage_data["prophages"])
        
        with open(output.json, 'w') as f:
            json.dump(prophage_data, f, indent=2)

