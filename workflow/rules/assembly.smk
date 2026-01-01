"""
============================================================================
Module 2: Assembly
============================================================================
Reconstruct genome(s) from reads

Assemblers:
    - SPAdes: Short-read assembly (careful mode)
    - Flye: Long-read assembly
    - Unicycler: Hybrid assembly

Steps:
    1. Select assembler based on read mode
    2. Run assembly
    3. Optional polishing (Medaka/Pilon)
    4. Filter contigs by length
    5. QUAST quality assessment
"""

# ============================================================================
# Helper Functions
# ============================================================================

def get_assembly_input_short_r1(wildcards):
    """Get clean short reads R1 for assembly"""
    return f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_R1.clean.fastq.gz"

def get_assembly_input_short_r2(wildcards):
    """Get clean short reads R2 for assembly"""
    return f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_R2.clean.fastq.gz"

def get_assembly_input_long(wildcards):
    """Get clean long reads for assembly"""
    return f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_long.clean.fastq.gz"


# ============================================================================
# SPAdes Assembly (Short Reads)
# ============================================================================

rule spades_assembly:
    """Assemble short reads with SPAdes"""
    input:
        r1 = get_assembly_input_short_r1,
        r2 = get_assembly_input_short_r2,
    output:
        contigs = f"{OUTPUT_DIR}/02_assembly/spades/contigs.fasta",
        scaffolds = f"{OUTPUT_DIR}/02_assembly/spades/scaffolds.fasta",
        graph = f"{OUTPUT_DIR}/02_assembly/spades/assembly_graph.fastg",
    params:
        outdir = f"{OUTPUT_DIR}/02_assembly/spades",
        mode = lambda w: get_assembly_params(config)["spades"]["mode"],
        kmers = lambda w: ",".join(map(str, get_assembly_params(config)["spades"]["kmer_sizes"])),
        memory = get_memory_gb(config),
    threads: workflow.cores
    log:
        f"{OUTPUT_DIR}/00_logs/spades.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        # Remove existing output dir if present (SPAdes requirement)
        rm -rf {params.outdir}
        
        spades.py \
            -1 {input.r1} \
            -2 {input.r2} \
            -o {params.outdir} \
            -t {threads} \
            -m {params.memory} \
            -k {params.kmers} \
            {"--careful" if params.mode == "careful" else ""} \
            2>&1 | tee {log}
        """


# ============================================================================
# Flye Assembly (Long Reads)
# ============================================================================

rule flye_assembly:
    """Assemble long reads with Flye"""
    input:
        get_assembly_input_long
    output:
        contigs = f"{OUTPUT_DIR}/02_assembly/flye/assembly.fasta",
        info = f"{OUTPUT_DIR}/02_assembly/flye/assembly_info.txt",
        graph = f"{OUTPUT_DIR}/02_assembly/flye/assembly_graph.gfa",
    params:
        outdir = f"{OUTPUT_DIR}/02_assembly/flye",
        read_type = lambda w: "--" + get_assembly_params(config)["flye"]["read_type"],
        genome_size = lambda w: get_assembly_params(config)["flye"]["genome_size"],
    threads: workflow.cores
    log:
        f"{OUTPUT_DIR}/00_logs/flye.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        flye \
            {params.read_type} {input} \
            --out-dir {params.outdir} \
            --threads {threads} \
            --genome-size {params.genome_size} \
            2>&1 | tee {log}
        """


# ============================================================================
# Unicycler Assembly (Hybrid)
# ============================================================================

rule unicycler_assembly:
    """Hybrid assembly with Unicycler"""
    input:
        r1 = get_assembly_input_short_r1,
        r2 = get_assembly_input_short_r2,
        long = get_assembly_input_long,
    output:
        contigs = f"{OUTPUT_DIR}/02_assembly/unicycler/assembly.fasta",
        graph = f"{OUTPUT_DIR}/02_assembly/unicycler/assembly.gfa",
        log_file = f"{OUTPUT_DIR}/02_assembly/unicycler/unicycler.log",
    params:
        outdir = f"{OUTPUT_DIR}/02_assembly/unicycler",
    threads: workflow.cores
    log:
        f"{OUTPUT_DIR}/00_logs/unicycler.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        unicycler \
            -1 {input.r1} \
            -2 {input.r2} \
            -l {input.long} \
            -o {params.outdir} \
            -t {threads} \
            2>&1 | tee {log}
        """


# ============================================================================
# Assembly Selection
# ============================================================================

def get_raw_assembly(wildcards):
    """Get the appropriate raw assembly based on read mode"""
    mode = detect_read_mode(config)
    assembler = select_assembler(config)
    
    if assembler == "spades" or mode == "short":
        return f"{OUTPUT_DIR}/02_assembly/spades/contigs.fasta"
    elif assembler == "flye" or mode == "long":
        return f"{OUTPUT_DIR}/02_assembly/flye/assembly.fasta"
    elif assembler == "unicycler" or mode == "hybrid":
        return f"{OUTPUT_DIR}/02_assembly/unicycler/assembly.fasta"
    else:
        return f"{OUTPUT_DIR}/02_assembly/spades/contigs.fasta"


rule select_assembly:
    """Select the appropriate assembly based on read mode"""
    input:
        get_raw_assembly
    output:
        f"{OUTPUT_DIR}/02_assembly/raw_assembly.fasta"
    shell:
        """
        cp {input} {output}
        """


# ============================================================================
# Polishing
# ============================================================================

rule medaka_polish:
    """Polish assembly with Medaka (Nanopore)"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/raw_assembly.fasta",
        reads = get_assembly_input_long,
    output:
        polished = f"{OUTPUT_DIR}/02_assembly/medaka/consensus.fasta",
    params:
        outdir = f"{OUTPUT_DIR}/02_assembly/medaka",
    threads: workflow.cores
    log:
        f"{OUTPUT_DIR}/00_logs/medaka.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        medaka_consensus \
            -i {input.reads} \
            -d {input.assembly} \
            -o {params.outdir} \
            -t {threads} \
            2>&1 | tee {log}
        """


rule pilon_polish:
    """Polish assembly with Pilon (Illumina)"""
    input:
        assembly = f"{OUTPUT_DIR}/02_assembly/raw_assembly.fasta",
        r1 = get_assembly_input_short_r1,
        r2 = get_assembly_input_short_r2,
    output:
        polished = f"{OUTPUT_DIR}/02_assembly/pilon/pilon.fasta",
        bam = f"{OUTPUT_DIR}/02_assembly/pilon/aligned.sorted.bam",
    params:
        outdir = f"{OUTPUT_DIR}/02_assembly/pilon",
    threads: workflow.cores
    log:
        f"{OUTPUT_DIR}/00_logs/pilon.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        
        # Index assembly
        bwa index {input.assembly}
        
        # Align reads
        bwa mem -t {threads} {input.assembly} {input.r1} {input.r2} | \
            samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        
        # Run Pilon
        pilon \
            --genome {input.assembly} \
            --frags {output.bam} \
            --outdir {params.outdir} \
            --output pilon \
            --threads {threads} \
            2>&1 | tee {log}
        """


def get_polished_assembly(wildcards):
    """Get polished or raw assembly based on config"""
    polish = get_assembly_params(config)["polish"]
    mode = detect_read_mode(config)
    
    if not polish:
        return f"{OUTPUT_DIR}/02_assembly/raw_assembly.fasta"
    
    if mode == "long":
        return f"{OUTPUT_DIR}/02_assembly/medaka/consensus.fasta"
    elif mode in ["short", "hybrid"]:
        return f"{OUTPUT_DIR}/02_assembly/pilon/pilon.fasta"
    else:
        return f"{OUTPUT_DIR}/02_assembly/raw_assembly.fasta"


rule select_polished:
    """Select polished or raw assembly"""
    input:
        get_polished_assembly
    output:
        f"{OUTPUT_DIR}/02_assembly/polished_assembly.fasta"
    shell:
        """
        cp {input} {output}
        """


# ============================================================================
# Filter Contigs
# ============================================================================

rule filter_contigs:
    """Filter contigs by minimum length"""
    input:
        f"{OUTPUT_DIR}/02_assembly/polished_assembly.fasta"
    output:
        f"{OUTPUT_DIR}/02_assembly/contigs.fasta"
    params:
        min_length = lambda w: get_assembly_params(config)["min_contig_length"],
    log:
        f"{OUTPUT_DIR}/00_logs/filter_contigs.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        seqkit seq -m {params.min_length} {input} > {output} 2> {log}
        """


# ============================================================================
# Assembly QC
# ============================================================================

rule quast:
    """Assess assembly quality with QUAST"""
    input:
        f"{OUTPUT_DIR}/02_assembly/contigs.fasta"
    output:
        report_html = f"{OUTPUT_DIR}/02_assembly/quast/report.html",
        report_tsv = f"{OUTPUT_DIR}/02_assembly/quast/report.tsv",
        report_txt = f"{OUTPUT_DIR}/02_assembly/quast/report.txt",
    params:
        outdir = f"{OUTPUT_DIR}/02_assembly/quast",
        min_contig = lambda w: get_assembly_params(config)["min_contig_length"],
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/quast.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        """
        quast.py \
            {input} \
            -o {params.outdir} \
            --min-contig {params.min_contig} \
            -t {threads} \
            2>&1 | tee {log}
        """


rule assembly_stats:
    """Generate assembly statistics JSON"""
    input:
        contigs = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        quast = f"{OUTPUT_DIR}/02_assembly/quast/report.tsv",
    output:
        f"{OUTPUT_DIR}/02_assembly/assembly_stats.json"
    run:
        import json
        from Bio import SeqIO
        
        # Parse assembly
        contigs = list(SeqIO.parse(input.contigs, "fasta"))
        lengths = [len(c) for c in contigs]
        total_length = sum(lengths)
        gc_content = sum(str(c.seq).upper().count('G') + str(c.seq).upper().count('C') for c in contigs) / total_length * 100 if total_length > 0 else 0
        
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
            "assembler": select_assembler(config),
            "num_contigs": len(contigs),
            "total_length": total_length,
            "largest_contig": max(lengths) if lengths else 0,
            "smallest_contig": min(lengths) if lengths else 0,
            "n50": n50,
            "gc_content": round(gc_content, 2),
            "mean_contig_length": round(total_length / len(contigs), 2) if contigs else 0,
        }
        
        with open(output[0], 'w') as f:
            json.dump(stats, f, indent=2)

