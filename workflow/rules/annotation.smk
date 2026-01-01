"""
============================================================================
Module 4: Annotation
============================================================================
Predict and annotate genomic features

Annotators:
    - Bakta: Bacterial annotation (comprehensive functional annotation)
    - Pharokka: Viral/phage annotation

Steps:
    1. Determine organism type
    2. Run appropriate annotator
    3. Extract proteins for downstream analysis
"""

# ============================================================================
# Helper Functions
# ============================================================================

def get_organism_type():
    """Get organism type from identification"""
    org_file = f"{OUTPUT_DIR}/03_identification/organism_type.txt"
    if os.path.exists(org_file):
        with open(org_file) as f:
            return f.read().strip()
    return "bacteria"  # Default


# ============================================================================
# Bakta Annotation (Bacteria)
# ============================================================================

rule bakta_annotate:
    """Annotate bacterial genome with Bakta"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        organism_type = f"{OUTPUT_DIR}/03_identification/organism_type.txt",
    output:
        gff = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3",
        gbk = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gbff",
        faa = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.faa",
        fna = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.fna",
        tsv = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.tsv",
        json = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.json",
        txt = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.txt",
    params:
        db = lambda w: get_database_path(config, "bakta"),
        outdir = f"{OUTPUT_DIR}/04_annotation/bakta",
        prefix = SAMPLE_ID,
        min_contig = lambda w: config.get("annotation", {}).get("bakta", {}).get("min_contig_length", 1),
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/bakta.log"
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        bakta \
            --db {params.db} \
            --output {params.outdir} \
            --prefix {params.prefix} \
            --min-contig-length {params.min_contig} \
            --threads {threads} \
            --force \
            {input.assembly} \
            2>&1 | tee {log}
        """


# ============================================================================
# Pharokka Annotation (Viruses/Phages)
# ============================================================================

rule pharokka_annotate:
    """Annotate viral/phage genome with Pharokka"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        organism_type = f"{OUTPUT_DIR}/03_identification/organism_type.txt",
    output:
        gff = f"{OUTPUT_DIR}/04_annotation/pharokka/{SAMPLE_ID}.gff",
        gbk = f"{OUTPUT_DIR}/04_annotation/pharokka/{SAMPLE_ID}.gbk",
        faa = f"{OUTPUT_DIR}/04_annotation/pharokka/{SAMPLE_ID}.faa",
        flag = f"{OUTPUT_DIR}/04_annotation/pharokka/pharokka_complete.flag",
    params:
        db = lambda w: get_database_path(config, "pharokka"),
        outdir = f"{OUTPUT_DIR}/04_annotation/pharokka",
        prefix = SAMPLE_ID,
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/pharokka.log"
    conda:
        "../envs/annotation.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        pharokka.py \
            --infile {input.assembly} \
            --outdir {params.outdir} \
            --database {params.db} \
            --prefix {params.prefix} \
            --threads {threads} \
            --force \
            2>&1 | tee {log}
        
        touch {output.flag}
        """


# ============================================================================
# Annotation Selection
# ============================================================================

def get_annotation_outputs(wildcards):
    """Get appropriate annotation outputs based on organism type"""
    org_file = f"{OUTPUT_DIR}/03_identification/organism_type.txt"
    
    # Default to bacteria if we can't determine
    org_type = "bacteria"
    if os.path.exists(org_file):
        with open(org_file) as f:
            org_type = f.read().strip()
    
    if org_type == "virus":
        return {
            "gff": f"{OUTPUT_DIR}/04_annotation/pharokka/{SAMPLE_ID}.gff",
            "proteins": f"{OUTPUT_DIR}/04_annotation/pharokka/{SAMPLE_ID}.faa",
        }
    else:
        return {
            "gff": f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3",
            "proteins": f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.faa",
        }


rule select_annotation:
    """Select appropriate annotation based on organism type"""
    input:
        organism_type = f"{OUTPUT_DIR}/03_identification/organism_type.txt",
        bakta_gff = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3",
        bakta_faa = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.faa",
    output:
        gff = f"{OUTPUT_DIR}/04_annotation/annotation.gff",
        proteins = f"{OUTPUT_DIR}/04_annotation/proteins.faa",
        flag = f"{OUTPUT_DIR}/04_annotation/annotation_complete.flag",
    run:
        import shutil
        
        # Read organism type
        with open(input.organism_type) as f:
            org_type = f.read().strip()
        
        # Copy appropriate files
        if org_type == "virus":
            pharokka_gff = f"{OUTPUT_DIR}/04_annotation/pharokka/{SAMPLE_ID}.gff"
            pharokka_faa = f"{OUTPUT_DIR}/04_annotation/pharokka/{SAMPLE_ID}.faa"
            if os.path.exists(pharokka_gff):
                shutil.copy(pharokka_gff, output.gff)
            if os.path.exists(pharokka_faa):
                shutil.copy(pharokka_faa, output.proteins)
        else:
            shutil.copy(input.bakta_gff, output.gff)
            shutil.copy(input.bakta_faa, output.proteins)
        
        # Touch flag
        with open(output.flag, 'w') as f:
            f.write(f"Annotation complete: {org_type}\n")


# ============================================================================
# Annotation Statistics
# ============================================================================

rule annotation_stats:
    """Generate annotation statistics"""
    input:
        gff = f"{OUTPUT_DIR}/04_annotation/annotation.gff",
        proteins = f"{OUTPUT_DIR}/04_annotation/proteins.faa",
        organism_type = f"{OUTPUT_DIR}/03_identification/organism_type.txt",
    output:
        f"{OUTPUT_DIR}/04_annotation/annotation_stats.json"
    run:
        import json
        from Bio import SeqIO
        
        # Read organism type
        with open(input.organism_type) as f:
            org_type = f.read().strip()
        
        # Count features from GFF
        feature_counts = {}
        with open(input.gff) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    feature_type = parts[2]
                    feature_counts[feature_type] = feature_counts.get(feature_type, 0) + 1
        
        # Count proteins
        proteins = list(SeqIO.parse(input.proteins, "fasta"))
        protein_lengths = [len(p) for p in proteins]
        
        stats = {
            "sample_id": SAMPLE_ID,
            "organism_type": org_type,
            "annotator": "pharokka" if org_type == "virus" else "bakta",
            "total_genes": feature_counts.get("gene", 0) or feature_counts.get("CDS", 0),
            "cds_count": feature_counts.get("CDS", 0),
            "trna_count": feature_counts.get("tRNA", 0),
            "rrna_count": feature_counts.get("rRNA", 0),
            "protein_count": len(proteins),
            "mean_protein_length": round(sum(protein_lengths) / len(protein_lengths), 1) if protein_lengths else 0,
            "feature_counts": feature_counts,
        }
        
        with open(output[0], 'w') as f:
            json.dump(stats, f, indent=2)

