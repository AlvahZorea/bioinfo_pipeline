"""
============================================================================
Module 1: Quality Control
============================================================================
Assess read quality, trim/filter reads, detect and remove contamination

Steps:
    1. FastQC on raw reads
    2. Trimming: fastp (Illumina) or Porechop+NanoFilt (Nanopore)
    3. Kraken2 contamination screening
    4. Remove contaminant reads (optional)
    5. FastQC on filtered reads
    6. MultiQC aggregate report
"""

# ============================================================================
# Helper Functions
# ============================================================================

def get_qc_input_r1(wildcards):
    """Get input R1 for QC"""
    r1, _ = get_short_reads(config)
    return r1

def get_qc_input_r2(wildcards):
    """Get input R2 for QC"""
    _, r2 = get_short_reads(config)
    return r2

def get_qc_input_long(wildcards):
    """Get input long reads for QC"""
    return get_long_reads(config)

# ============================================================================
# Illumina QC Rules
# ============================================================================

rule fastqc_raw_illumina:
    """Run FastQC on raw Illumina reads"""
    input:
        r1 = get_qc_input_r1,
        r2 = get_qc_input_r2,
    output:
        html_r1 = f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_R1_fastqc.html",
        html_r2 = f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_R2_fastqc.html",
        zip_r1 = f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_R1_fastqc.zip",
        zip_r2 = f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_R2_fastqc.zip",
    params:
        outdir = f"{OUTPUT_DIR}/01_qc/fastqc_raw"
    threads: 2
    log:
        f"{OUTPUT_DIR}/00_logs/fastqc_raw.log"
    conda:
        "../envs/qc.yaml"
    params:
        outdir = f"{OUTPUT_DIR}/01_qc/fastqc_raw",
        sample_id = SAMPLE_ID,
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2>&1 | tee {log}
        
        # Rename outputs to standard names
        cd {params.outdir}
        SAMPLE_ID="{params.sample_id}"
        for f in *_fastqc.html *_fastqc.zip; do
            if [[ "$f" == *"_1"* ]] || [[ "$f" == *"_R1"* ]]; then
                ext="${{f##*.}}"
                newname="${{SAMPLE_ID}}_R1_fastqc.${{ext}}"
                [ "$f" != "$newname" ] && mv "$f" "$newname" 2>/dev/null || true
            elif [[ "$f" == *"_2"* ]] || [[ "$f" == *"_R2"* ]]; then
                ext="${{f##*.}}"
                newname="${{SAMPLE_ID}}_R2_fastqc.${{ext}}"
                [ "$f" != "$newname" ] && mv "$f" "$newname" 2>/dev/null || true
            fi
        done
        """


rule fastp_trim:
    """Trim and filter Illumina reads with fastp"""
    input:
        r1 = get_qc_input_r1,
        r2 = get_qc_input_r2,
    output:
        r1 = f"{OUTPUT_DIR}/01_qc/trimmed/{SAMPLE_ID}_R1.trimmed.fastq.gz",
        r2 = f"{OUTPUT_DIR}/01_qc/trimmed/{SAMPLE_ID}_R2.trimmed.fastq.gz",
        html = f"{OUTPUT_DIR}/01_qc/trimmed/{SAMPLE_ID}_fastp.html",
        json = f"{OUTPUT_DIR}/01_qc/trimmed/{SAMPLE_ID}_fastp.json",
    params:
        min_quality = lambda w: get_qc_params(config)["min_quality"],
        min_length = lambda w: get_qc_params(config)["min_length"],
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/fastp.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.r1})
        fastp \
            -i {input.r1} \
            -I {input.r2} \
            -o {output.r1} \
            -O {output.r2} \
            --html {output.html} \
            --json {output.json} \
            --qualified_quality_phred {params.min_quality} \
            --length_required {params.min_length} \
            --detect_adapter_for_pe \
            --correction \
            --thread {threads} \
            2>&1 | tee {log}
        """


rule kraken2_contamination_illumina:
    """Screen for contamination using Kraken2"""
    input:
        r1 = rules.fastp_trim.output.r1,
        r2 = rules.fastp_trim.output.r2,
    output:
        report = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_kraken2_report.txt",
        output = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_kraken2_output.txt",
        classified_r1 = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_classified_R1.fastq.gz",
        classified_r2 = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_classified_R2.fastq.gz",
        unclassified_r1 = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_unclassified_R1.fastq.gz",
        unclassified_r2 = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_unclassified_R2.fastq.gz",
    params:
        db = lambda w: get_qc_params(config)["contaminant_db"] or get_database_path(config, "kraken2", "standard"),
        classified_prefix = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_classified#.fastq",
        unclassified_prefix = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_unclassified#.fastq",
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/kraken2_contamination.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.report})
        
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --report {output.report} \
            --output {output.output} \
            --classified-out {params.classified_prefix} \
            --unclassified-out {params.unclassified_prefix} \
            {input.r1} {input.r2} \
            2>&1 | tee {log}
        
        # Compress outputs
        gzip -f {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_classified_1.fastq && \
            mv {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_classified_1.fastq.gz {output.classified_r1} || true
        gzip -f {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_classified_2.fastq && \
            mv {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_classified_2.fastq.gz {output.classified_r2} || true
        gzip -f {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_unclassified_1.fastq && \
            mv {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_unclassified_1.fastq.gz {output.unclassified_r1} || true
        gzip -f {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_unclassified_2.fastq && \
            mv {OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_unclassified_2.fastq.gz {output.unclassified_r2} || true
        """


rule filter_contamination_illumina:
    """
    Produce clean reads - either filtered (contaminants removed) or just trimmed
    based on config setting.
    """
    input:
        trimmed_r1 = rules.fastp_trim.output.r1,
        trimmed_r2 = rules.fastp_trim.output.r2,
        unclassified_r1 = rules.kraken2_contamination_illumina.output.unclassified_r1,
        unclassified_r2 = rules.kraken2_contamination_illumina.output.unclassified_r2,
        kraken_report = rules.kraken2_contamination_illumina.output.report,
    output:
        r1 = f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_R1.clean.fastq.gz",
        r2 = f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_R2.clean.fastq.gz",
        contaminant_report = f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_contamination_summary.txt",
    params:
        remove_contaminants = lambda w: get_qc_params(config)["remove_contaminants"],
    log:
        f"{OUTPUT_DIR}/00_logs/filter_contamination.log"
    run:
        import shutil
        
        # Parse Kraken report to summarize contamination
        contaminants = []
        with open(input.kraken_report) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    pct, num_clade, num_direct, rank, taxid, name = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
                    if float(pct) > 0.1 and rank in ['S', 'G']:  # Species or Genus with >0.1%
                        contaminants.append(f"{name.strip()}: {pct}%")
        
        # Write contamination summary
        with open(output.contaminant_report, 'w') as f:
            f.write(f"# Contamination Summary for {SAMPLE_ID}\n")
            f.write(f"# Contaminants removed: {params.remove_contaminants}\n\n")
            if contaminants:
                f.write("Detected organisms (>0.1%):\n")
                for c in contaminants[:20]:  # Top 20
                    f.write(f"  {c}\n")
            else:
                f.write("No significant contamination detected.\n")
        
        # Copy appropriate files
        if params.remove_contaminants:
            # Use unclassified (non-contaminant) reads
            shutil.copy(input.unclassified_r1, output.r1)
            shutil.copy(input.unclassified_r2, output.r2)
        else:
            # Use all trimmed reads
            shutil.copy(input.trimmed_r1, output.r1)
            shutil.copy(input.trimmed_r2, output.r2)


rule fastqc_clean_illumina:
    """Run FastQC on clean reads"""
    input:
        r1 = rules.filter_contamination_illumina.output.r1,
        r2 = rules.filter_contamination_illumina.output.r2,
    output:
        html_r1 = f"{OUTPUT_DIR}/01_qc/fastqc_clean/{SAMPLE_ID}_R1.clean_fastqc.html",
        html_r2 = f"{OUTPUT_DIR}/01_qc/fastqc_clean/{SAMPLE_ID}_R2.clean_fastqc.html",
    params:
        outdir = f"{OUTPUT_DIR}/01_qc/fastqc_clean"
    threads: 2
    log:
        f"{OUTPUT_DIR}/00_logs/fastqc_clean.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input.r1} {input.r2} 2>&1 | tee {log}
        """


# ============================================================================
# Nanopore QC Rules
# ============================================================================

rule fastqc_raw_nanopore:
    """Run FastQC on raw Nanopore reads"""
    input:
        get_qc_input_long
    output:
        html = f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_long_fastqc.html",
        zip = f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_long_fastqc.zip",
    params:
        outdir = f"{OUTPUT_DIR}/01_qc/fastqc_raw"
    threads: 2
    log:
        f"{OUTPUT_DIR}/00_logs/fastqc_raw_long.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        fastqc -t {threads} -o {params.outdir} {input} 2>&1 | tee {log}
        """


rule porechop_trim:
    """Remove adapters from Nanopore reads with Porechop"""
    input:
        get_qc_input_long
    output:
        fastq = f"{OUTPUT_DIR}/01_qc/trimmed/{SAMPLE_ID}_long.porechop.fastq.gz",
    threads: 4
    log:
        f"{OUTPUT_DIR}/00_logs/porechop.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.fastq})
        porechop -i {input} -o {output.fastq} -t {threads} 2>&1 | tee {log}
        """


rule nanofilt_filter:
    """Quality filter Nanopore reads with NanoFilt"""
    input:
        rules.porechop_trim.output.fastq
    output:
        fastq = f"{OUTPUT_DIR}/01_qc/trimmed/{SAMPLE_ID}_long.filtered.fastq.gz",
    params:
        min_length = lambda w: get_qc_params(config)["nanopore"]["min_length"],
        min_quality = lambda w: get_qc_params(config)["nanopore"]["min_quality"],
    log:
        f"{OUTPUT_DIR}/00_logs/nanofilt.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        gunzip -c {input} | \
        NanoFilt -q {params.min_quality} -l {params.min_length} | \
        gzip > {output.fastq} 2>&1 | tee {log}
        """


rule kraken2_contamination_nanopore:
    """Screen Nanopore reads for contamination using Kraken2"""
    input:
        rules.nanofilt_filter.output.fastq
    output:
        report = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_long_kraken2_report.txt",
        output = f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_long_kraken2_output.txt",
    params:
        db = lambda w: get_qc_params(config)["contaminant_db"] or get_database_path(config, "kraken2", "standard"),
    threads: 8
    log:
        f"{OUTPUT_DIR}/00_logs/kraken2_contamination_long.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        mkdir -p $(dirname {output.report})
        
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --gzip-compressed \
            --report {output.report} \
            --output {output.output} \
            {input} \
            2>&1 | tee {log}
        """


rule filter_contamination_nanopore:
    """Produce clean Nanopore reads"""
    input:
        filtered = rules.nanofilt_filter.output.fastq,
        kraken_report = rules.kraken2_contamination_nanopore.output.report,
    output:
        fastq = f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_long.clean.fastq.gz",
        contaminant_report = f"{OUTPUT_DIR}/01_qc/clean/{SAMPLE_ID}_long_contamination_summary.txt",
    log:
        f"{OUTPUT_DIR}/00_logs/filter_contamination_long.log"
    run:
        import shutil
        
        # Parse Kraken report to summarize contamination
        contaminants = []
        with open(input.kraken_report) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 6:
                    pct, num_clade, num_direct, rank, taxid, name = parts[0], parts[1], parts[2], parts[3], parts[4], parts[5]
                    if float(pct) > 0.1 and rank in ['S', 'G']:
                        contaminants.append(f"{name.strip()}: {pct}%")
        
        # Write contamination summary
        with open(output.contaminant_report, 'w') as f:
            f.write(f"# Contamination Summary for {SAMPLE_ID} (long reads)\n\n")
            if contaminants:
                f.write("Detected organisms (>0.1%):\n")
                for c in contaminants[:20]:
                    f.write(f"  {c}\n")
            else:
                f.write("No significant contamination detected.\n")
        
        # Copy filtered reads (TODO: implement actual filtering if needed)
        shutil.copy(input.filtered, output.fastq)


# ============================================================================
# Aggregate QC
# ============================================================================

def get_multiqc_inputs(wildcards):
    """Get all QC outputs for MultiQC based on read mode"""
    inputs = []
    mode = detect_read_mode(config)
    
    if mode in ["short", "hybrid"]:
        inputs.extend([
            f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_R1_fastqc.zip",
            f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_R2_fastqc.zip",
            f"{OUTPUT_DIR}/01_qc/trimmed/{SAMPLE_ID}_fastp.json",
            f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_kraken2_report.txt",
            f"{OUTPUT_DIR}/01_qc/fastqc_clean/{SAMPLE_ID}_R1.clean_fastqc.html",
        ])
    
    if mode in ["long", "hybrid"]:
        inputs.extend([
            f"{OUTPUT_DIR}/01_qc/fastqc_raw/{SAMPLE_ID}_long_fastqc.zip",
            f"{OUTPUT_DIR}/01_qc/contamination/{SAMPLE_ID}_long_kraken2_report.txt",
        ])
    
    return inputs


rule multiqc:
    """Aggregate all QC reports with MultiQC"""
    input:
        get_multiqc_inputs
    output:
        html = f"{OUTPUT_DIR}/01_qc/multiqc_report.html",
        data = directory(f"{OUTPUT_DIR}/01_qc/multiqc_data"),
    params:
        qc_dir = f"{OUTPUT_DIR}/01_qc"
    log:
        f"{OUTPUT_DIR}/00_logs/multiqc.log"
    conda:
        "../envs/qc.yaml"
    shell:
        """
        multiqc \
            --force \
            --outdir {params.qc_dir} \
            --filename multiqc_report.html \
            {params.qc_dir} \
            2>&1 | tee {log}
        """

