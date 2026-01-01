"""
============================================================================
Module 8: Report Generation
============================================================================
Generate comprehensive HTML report with visualizations

Report sections:
    1. Executive Summary
    2. QC Results
    3. Assembly Statistics
    4. Taxonomic Identification
    5. Annotation Summary
    6. Genome Visualization (circular + linear alignment)
    7. Anomaly Detection (foreign genes, structural variants)
    8. AMR Analysis
    9. Virulence Factors
    10. Mobile Elements
    11. Prophages
"""

# ============================================================================
# Aggregate Results
# ============================================================================

rule aggregate_results:
    """Aggregate all analysis results into a single JSON"""
    input:
        # Assembly
        assembly_stats = f"{OUTPUT_DIR}/02_assembly/assembly_stats.json",
        contigs = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        # Identification
        taxonomy = f"{OUTPUT_DIR}/03_identification/taxonomy.json",
        organism_type = f"{OUTPUT_DIR}/03_identification/organism_type.txt",
        # Annotation
        annotation_stats = f"{OUTPUT_DIR}/04_annotation/annotation_stats.json",
        annotation_flag = f"{OUTPUT_DIR}/04_annotation/annotation_complete.flag",
        # Specialized
        amr = f"{OUTPUT_DIR}/05_specialized/amr/amr_summary.json",
        virulence = f"{OUTPUT_DIR}/05_specialized/virulence/virulence_summary.json",
        mge = f"{OUTPUT_DIR}/05_specialized/mge/mge_summary.json",
        prophage = f"{OUTPUT_DIR}/05_specialized/prophages/prophage_summary.json",
    output:
        aggregated = f"{OUTPUT_DIR}/08_report/all_results.json",
    run:
        import json
        import datetime
        from Bio import SeqIO
        
        results = {
            "sample_id": SAMPLE_ID,
            "report_generated": datetime.datetime.now().isoformat(),
            "pipeline_version": "1.0.0",
        }
        
        # Load all JSON results
        json_inputs = {
            "assembly": input.assembly_stats,
            "taxonomy": input.taxonomy,
            "annotation": input.annotation_stats,
            "amr": input.amr,
            "virulence": input.virulence,
            "mge": input.mge,
            "prophages": input.prophage,
        }
        
        # Add anomaly data if available
        anomaly_path = f"{OUTPUT_DIR}/06_anomaly/anomaly_summary.json"
        if os.path.exists(anomaly_path):
            json_inputs["anomaly"] = anomaly_path
        
        # Add synteny data if available
        synteny_path = f"{OUTPUT_DIR}/06_anomaly/synteny/synteny_report.json"
        if os.path.exists(synteny_path):
            json_inputs["synteny"] = synteny_path
        
        # Add intergenic data if available
        intergenic_path = f"{OUTPUT_DIR}/06_anomaly/intergenic/intergenic_taxonomy.json"
        if os.path.exists(intergenic_path):
            json_inputs["intergenic"] = intergenic_path
        
        # Add integrated elements data if available
        integrated_path = f"{OUTPUT_DIR}/06_anomaly/integrated_elements/integrated_elements.json"
        if os.path.exists(integrated_path):
            json_inputs["integrated_elements"] = integrated_path
        
        for key, filepath in json_inputs.items():
            if os.path.exists(filepath):
                with open(filepath) as f:
                    results[key] = json.load(f)
        
        # Read organism type
        with open(input.organism_type) as f:
            results["organism_type"] = f.read().strip()
        
        # Load contig info for visualization
        contigs_data = []
        for record in SeqIO.parse(input.contigs, "fasta"):
            seq = str(record.seq).upper()
            gc = (seq.count('G') + seq.count('C')) / len(seq) * 100 if seq else 0
            contigs_data.append({
                "id": record.id,
                "length": len(record),
                "gc": round(gc, 2)
            })
        results["contigs"] = contigs_data
        
        # Write aggregated results
        with open(output.aggregated, 'w') as f:
            json.dump(results, f, indent=2)


# ============================================================================
# Report Generation
# ============================================================================

rule generate_report:
    """Generate HTML report with genome visualization"""
    input:
        aggregated = rules.aggregate_results.output.aggregated,
        contigs = f"{OUTPUT_DIR}/02_assembly/contigs.fasta",
        quast_report = f"{OUTPUT_DIR}/02_assembly/quast/report.html",
        multiqc_report = f"{OUTPUT_DIR}/01_qc/multiqc_report.html",
    output:
        html = f"{OUTPUT_DIR}/08_report/genome_report.html",
    params:
        self_contained = lambda w: config.get("report", {}).get("self_contained", True),
        title = lambda w: config.get("report", {}).get("title") or f"Genome Analysis Report: {SAMPLE_ID}",
    log:
        f"{OUTPUT_DIR}/00_logs/generate_report.log"
    run:
        import json
        from Bio import SeqIO
        
        # Load aggregated results
        with open(input.aggregated) as f:
            results = json.load(f)
        
        # Load annotation features if available (check file exists)
        features = []
        gff_path = f"{OUTPUT_DIR}/04_annotation/bakta/{SAMPLE_ID}.gff3"
        if os.path.exists(gff_path):
            features = parse_gff_features(gff_path)
        results["features"] = features
        
        # Load alignment data if available (check file exists)
        alignment_regions = []
        coords_path = f"{OUTPUT_DIR}/06_anomaly/alignment/nucmer.coords"
        if os.path.exists(coords_path):
            alignment_regions = parse_nucmer_coords(coords_path)
        results["alignment_regions"] = alignment_regions
        
        # Generate HTML report
        html_content = generate_html_report(results, params.title, input.contigs)
        
        with open(output.html, 'w') as f:
            f.write(html_content)


def parse_gff_features(gff_file):
    """Parse GFF3 annotation file to extract features for visualization"""
    features = []
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 9:
                ftype = parts[2]
                if ftype in ['CDS', 'gene', 'tRNA', 'rRNA']:
                    attrs = dict(x.split('=') for x in parts[8].split(';') if '=' in x)
                    features.append({
                        "contig": parts[0],
                        "type": ftype,
                        "start": int(parts[3]),
                        "end": int(parts[4]),
                        "strand": parts[6],
                        "name": attrs.get('Name', attrs.get('locus_tag', '')),
                        "product": attrs.get('product', '')
                    })
    return features


def parse_nucmer_coords(coords_file):
    """Parse nucmer coords file for alignment visualization"""
    regions = []
    with open(coords_file) as f:
        for line in f:
            if line.startswith('=') or line.startswith('[') or '|' in line:
                continue
            parts = line.strip().split()
            if len(parts) >= 12:
                try:
                    regions.append({
                        "ref_start": int(parts[0]),
                        "ref_end": int(parts[1]),
                        "qry_start": int(parts[3]),
                        "qry_end": int(parts[4]),
                        "identity": float(parts[9]),
                        "contig": parts[11]
                    })
                except (ValueError, IndexError):
                    pass
    return regions


def generate_html_report(results, title, contigs_file):
    """Generate self-contained HTML report with genome visualization"""
    from Bio import SeqIO
    import json
    
    # Load contig data for visualization
    contigs = list(SeqIO.parse(contigs_file, "fasta"))
    contig_lengths = [len(c) for c in contigs]
    contigs_data = results.get("contigs", [])
    features = results.get("features", [])
    alignment_regions = results.get("alignment_regions", [])
    anomaly = results.get("anomaly", {})
    
    # Check for anomalies
    has_anomaly = anomaly.get("anomaly_detected", False)
    foreign_genes = anomaly.get("gene_taxonomy", {}).get("foreign_genes", [])
    novel_genes = anomaly.get("gene_taxonomy", {}).get("novel_genes", [])
    foreign_clusters = anomaly.get("foreign_gene_clusters", [])
    synteny = results.get("synteny", {})
    intergenic = results.get("intergenic", {})
    integrated_elements = results.get("integrated_elements", {})
    minority_species = results.get("taxonomy", {}).get("kraken2", {}).get("minority_species", [])
    deleted_genes = anomaly.get("alignment", {}).get("deleted_genes", [])
    
    # Build report HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --success-color: #27ae60;
            --warning-color: #f39c12;
            --danger-color: #e74c3c;
            --anomaly-color: #9b59b6;
            --light-bg: #f8f9fa;
            --border-color: #dee2e6;
        }}
        
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: var(--light-bg);
        }}
        
        .container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
        
        header {{
            background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
            color: white;
            padding: 30px;
            margin-bottom: 30px;
            border-radius: 8px;
        }}
        
        header h1 {{ font-size: 2rem; margin-bottom: 10px; }}
        header .meta {{ opacity: 0.9; font-size: 0.9rem; }}
        
        .section {{
            background: white;
            border-radius: 8px;
            padding: 25px;
            margin-bottom: 25px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        
        .section h2 {{
            color: var(--primary-color);
            border-bottom: 2px solid var(--secondary-color);
            padding-bottom: 10px;
            margin-bottom: 20px;
        }}
        
        .section h3 {{ color: var(--primary-color); margin: 20px 0 10px; }}
        
        .section.anomaly-alert {{
            border: 3px solid var(--danger-color);
            background: linear-gradient(135deg, #fff5f5 0%, #ffe6e6 100%);
        }}
        
        .section.anomaly-alert h2 {{
            color: var(--danger-color);
            border-bottom-color: var(--danger-color);
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-bottom: 20px;
        }}
        
        .stat-card {{
            background: var(--light-bg);
            padding: 15px;
            border-radius: 8px;
            text-align: center;
            border-left: 4px solid var(--secondary-color);
        }}
        
        .stat-card.anomaly {{ border-left-color: var(--danger-color); background: #ffe6e6; }}
        
        .stat-card .value {{ font-size: 1.8rem; font-weight: bold; color: var(--primary-color); }}
        .stat-card .label {{ font-size: 0.85rem; color: #666; text-transform: uppercase; }}
        
        table {{ width: 100%; border-collapse: collapse; margin: 15px 0; }}
        th, td {{ padding: 12px; text-align: left; border-bottom: 1px solid var(--border-color); }}
        th {{ background: var(--light-bg); font-weight: 600; color: var(--primary-color); }}
        tr:hover {{ background: #f5f5f5; }}
        tr.foreign-gene {{ background: #ffe6e6; }}
        tr.foreign-gene:hover {{ background: #ffcccc; }}
        
        .badge {{
            display: inline-block;
            padding: 4px 10px;
            border-radius: 20px;
            font-size: 0.8rem;
            font-weight: 500;
        }}
        
        .badge-success {{ background: #d4edda; color: #155724; }}
        .badge-warning {{ background: #fff3cd; color: #856404; }}
        .badge-danger {{ background: #f8d7da; color: #721c24; }}
        .badge-info {{ background: #d1ecf1; color: #0c5460; }}
        .badge-anomaly {{ background: #e8daef; color: #6c3483; }}
        
        .plot-container {{ margin: 20px 0; min-height: 400px; }}
        .genome-viz-container {{ margin: 20px 0; min-height: 500px; }}
        
        .summary-box {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 20px;
        }}
        
        .summary-box.alert-box {{
            background: linear-gradient(135deg, #e74c3c 0%, #c0392b 100%);
        }}
        
        .summary-box h3 {{ color: white; margin-bottom: 15px; }}
        .summary-item {{ display: flex; justify-content: space-between; padding: 8px 0; border-bottom: 1px solid rgba(255,255,255,0.2); }}
        .summary-item:last-child {{ border-bottom: none; }}
        
        .alert {{ padding: 15px; border-radius: 8px; margin: 15px 0; }}
        .alert-info {{ background: #e7f3ff; border-left: 4px solid var(--secondary-color); }}
        .alert-warning {{ background: #fff8e6; border-left: 4px solid var(--warning-color); }}
        .alert-danger {{ background: #ffe6e6; border-left: 4px solid var(--danger-color); }}
        
        .legend {{ display: flex; flex-wrap: wrap; gap: 15px; margin: 15px 0; padding: 15px; background: #f8f9fa; border-radius: 8px; }}
        .legend-item {{ display: flex; align-items: center; gap: 8px; }}
        .legend-color {{ width: 20px; height: 20px; border-radius: 4px; }}
        
        footer {{ text-align: center; padding: 20px; color: #666; font-size: 0.9rem; }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>{title}</h1>
            <div class="meta">
                <strong>Sample ID:</strong> {results.get('sample_id', 'N/A')} | 
                <strong>Organism Type:</strong> {results.get('organism_type', 'Unknown')} | 
                <strong>Generated:</strong> {results.get('report_generated', 'N/A')[:19]}
                {"| <span style='color: #ffcccc;'>⚠️ ANOMALIES DETECTED</span>" if has_anomaly else ""}
            </div>
        </header>
        
        <!-- Anomaly Alert (if detected) -->
        {generate_anomaly_alert_section(results) if has_anomaly else ""}
        
        <!-- Executive Summary -->
        <div class="section">
            <h2>Executive Summary</h2>
            <div class="summary-box{' alert-box' if has_anomaly else ''}">
                <h3>Key Findings</h3>
                {generate_executive_summary(results)}
            </div>
        </div>
        
        <!-- Genome Visualization -->
        <div class="section">
            <h2>Genome Visualization</h2>
            <p>Interactive visualization of genome structure, features, and anomalies.</p>
            
            <h3>Circular Genome Map</h3>
            <div id="circular-genome-plot" class="genome-viz-container"></div>
            
            <div class="legend">
                <div class="legend-item"><div class="legend-color" style="background: #3498db;"></div> CDS</div>
                <div class="legend-item"><div class="legend-color" style="background: #27ae60;"></div> tRNA</div>
                <div class="legend-item"><div class="legend-color" style="background: #f39c12;"></div> rRNA</div>
                <div class="legend-item"><div class="legend-color" style="background: #e74c3c;"></div> Foreign Gene</div>
                <div class="legend-item"><div class="legend-color" style="background: #9b59b6;"></div> Anomalous Region</div>
            </div>
            
            <h3>Linear Contig View</h3>
            <div id="linear-contig-plot" class="genome-viz-container"></div>
        </div>
        
        <!-- Reference Alignment (if available) -->
        {generate_alignment_section(results)}
        
        <!-- Assembly Statistics -->
        <div class="section">
            <h2>Assembly Statistics</h2>
            {generate_assembly_section(results)}
            <div id="contig-length-plot" class="plot-container"></div>
        </div>
        
        <!-- Taxonomy -->
        <div class="section">
            <h2>Taxonomic Identification</h2>
            {generate_taxonomy_section(results)}
            <div id="taxonomy-pie-plot" class="plot-container"></div>
        </div>
        
        <!-- Anomaly Detection Details -->
        {generate_anomaly_details_section(results)}
        
        <!-- Annotation -->
        <div class="section">
            <h2>Genome Annotation</h2>
            {generate_annotation_section(results)}
        </div>
        
        <!-- AMR -->
        <div class="section">
            <h2>Antimicrobial Resistance</h2>
            {generate_amr_section(results)}
        </div>
        
        <!-- Virulence -->
        <div class="section">
            <h2>Virulence Factors</h2>
            {generate_virulence_section(results)}
        </div>
        
        <!-- Mobile Elements -->
        <div class="section">
            <h2>Mobile Genetic Elements</h2>
            {generate_mge_section(results)}
        </div>
        
        <!-- Prophages -->
        <div class="section">
            <h2>Prophage Regions</h2>
            {generate_prophage_section(results)}
        </div>
        
        <footer>
            Genome Analysis Pipeline v{results.get('pipeline_version', '1.0.0')} | Report generated: {results.get('report_generated', 'N/A')[:19]}
        </footer>
    </div>
    
    <script>
        // Data
        var contigsData = {json.dumps(contigs_data)};
        var features = {json.dumps(features[:500])};  // Limit for performance
        var alignmentRegions = {json.dumps(alignment_regions)};
        var foreignGenes = {json.dumps(foreign_genes)};
        var novelGenes = {json.dumps(novel_genes)};
        var foreignClusters = {json.dumps(foreign_clusters)};
        var deletedGenes = {json.dumps(deleted_genes)};
        var taxonomyDist = {json.dumps(anomaly.get('gene_taxonomy', {}).get('taxonomy_distribution', {}))};
        var syntenyData = {json.dumps(synteny)};
        var intergenicData = {json.dumps(intergenic)};
        var integratedElements = {json.dumps(integrated_elements.get('candidates', []) if isinstance(integrated_elements, dict) else [])};
        
        // Circular Genome Map
        function drawCircularGenome() {{
            var totalLength = contigsData.reduce((sum, c) => sum + c.length, 0);
            var traces = [];
            
            // Outer ring - contigs
            var theta = [];
            var r = [];
            var colors = [];
            var texts = [];
            
            var cumLength = 0;
            contigsData.forEach((contig, i) => {{
                var startAngle = (cumLength / totalLength) * 360;
                var endAngle = ((cumLength + contig.length) / totalLength) * 360;
                
                for (var a = startAngle; a <= endAngle; a += 1) {{
                    theta.push(a);
                    r.push(1);
                    colors.push(contig.gc);
                    texts.push(contig.id + '<br>Length: ' + contig.length.toLocaleString() + ' bp<br>GC: ' + contig.gc + '%');
                }}
                cumLength += contig.length;
            }});
            
            traces.push({{
                type: 'scatterpolar',
                mode: 'lines',
                r: r,
                theta: theta,
                line: {{ width: 15 }},
                marker: {{ color: colors, colorscale: 'RdYlBu', reversescale: true, showscale: true, colorbar: {{ title: 'GC%' }} }},
                text: texts,
                hoverinfo: 'text',
                name: 'Contigs'
            }});
            
            // Inner rings - features (simplified)
            if (features.length > 0) {{
                var cdsTheta = [];
                var cdsR = [];
                var cdsTexts = [];
                var foreignTheta = [];
                var foreignR = [];
                var foreignTexts = [];
                
                features.forEach(f => {{
                    var contigIdx = contigsData.findIndex(c => c.id === f.contig);
                    if (contigIdx >= 0) {{
                        var contigStart = contigsData.slice(0, contigIdx).reduce((sum, c) => sum + c.length, 0);
                        var angle = ((contigStart + f.start) / totalLength) * 360;
                        
                        var isForeign = foreignGenes.some(fg => fg.gene_id === f.name);
                        
                        if (isForeign) {{
                            foreignTheta.push(angle);
                            foreignR.push(0.75);
                            foreignTexts.push(f.name + ' (FOREIGN)<br>' + f.product);
                        }} else if (f.type === 'CDS') {{
                            cdsTheta.push(angle);
                            cdsR.push(0.85);
                            cdsTexts.push(f.name + '<br>' + f.product);
                        }}
                    }}
                }});
                
                if (cdsTheta.length > 0) {{
                    traces.push({{
                        type: 'scatterpolar',
                        mode: 'markers',
                        r: cdsR,
                        theta: cdsTheta,
                        marker: {{ size: 3, color: '#3498db' }},
                        text: cdsTexts,
                        hoverinfo: 'text',
                        name: 'CDS'
                    }});
                }}
                
                if (foreignTheta.length > 0) {{
                    traces.push({{
                        type: 'scatterpolar',
                        mode: 'markers',
                        r: foreignR,
                        theta: foreignTheta,
                        marker: {{ size: 8, color: '#e74c3c', symbol: 'diamond' }},
                        text: foreignTexts,
                        hoverinfo: 'text',
                        name: 'Foreign Genes'
                    }});
                }}
            }}
            
            var layout = {{
                polar: {{
                    radialaxis: {{ visible: false, range: [0, 1.2] }},
                    angularaxis: {{ direction: 'clockwise' }}
                }},
                showlegend: true,
                title: 'Circular Genome Map (GC content colored)',
                height: 600
            }};
            
            Plotly.newPlot('circular-genome-plot', traces, layout);
        }}
        
        // Linear Contig View
        function drawLinearContigs() {{
            var traces = [];
            var shapes = [];
            var annotations = [];
            
            contigsData.forEach((contig, i) => {{
                // Contig bar
                shapes.push({{
                    type: 'rect',
                    x0: 0,
                    x1: contig.length,
                    y0: i - 0.3,
                    y1: i + 0.3,
                    fillcolor: '#3498db',
                    opacity: 0.6,
                    line: {{ width: 1, color: '#2c3e50' }}
                }});
                
                annotations.push({{
                    x: -5000,
                    y: i,
                    text: contig.id,
                    showarrow: false,
                    xanchor: 'right'
                }});
            }});
            
            // Foreign genes markers
            if (foreignGenes.length > 0 && features.length > 0) {{
                var foreignX = [];
                var foreignY = [];
                var foreignText = [];
                
                foreignGenes.forEach(fg => {{
                    var feature = features.find(f => f.name === fg.gene_id);
                    if (feature) {{
                        var contigIdx = contigsData.findIndex(c => c.id === feature.contig);
                        if (contigIdx >= 0) {{
                            foreignX.push((feature.start + feature.end) / 2);
                            foreignY.push(contigIdx);
                            foreignText.push(fg.gene_id + '<br>Origin: ' + fg.organism + '<br>Identity: ' + fg.identity.toFixed(1) + '%');
                        }}
                    }}
                }});
                
                if (foreignX.length > 0) {{
                    traces.push({{
                        type: 'scatter',
                        mode: 'markers',
                        x: foreignX,
                        y: foreignY,
                        marker: {{ size: 12, color: '#e74c3c', symbol: 'triangle-up' }},
                        text: foreignText,
                        hoverinfo: 'text',
                        name: 'Foreign Genes'
                    }});
                }}
            }}
            
            var layout = {{
                shapes: shapes,
                annotations: annotations,
                xaxis: {{ title: 'Position (bp)', rangemode: 'tozero' }},
                yaxis: {{ 
                    tickvals: contigsData.map((c, i) => i),
                    ticktext: contigsData.map(c => ''),
                    range: [-1, contigsData.length]
                }},
                title: 'Linear Contig View with Anomalies',
                height: Math.max(300, contigsData.length * 50 + 100),
                showlegend: true
            }};
            
            Plotly.newPlot('linear-contig-plot', traces, layout);
        }}
        
        // Contig length histogram
        var contigLengths = {json.dumps(contig_lengths)};
        Plotly.newPlot('contig-length-plot', [{{
            x: contigLengths,
            type: 'histogram',
            marker: {{ color: '#3498db' }},
            nbinsx: 50
        }}], {{
            title: 'Contig Length Distribution',
            xaxis: {{ title: 'Contig Length (bp)' }},
            yaxis: {{ title: 'Count' }},
            bargap: 0.05
        }});
        
        // Taxonomy pie chart
        if (Object.keys(taxonomyDist).length > 0) {{
            var labels = Object.keys(taxonomyDist).slice(0, 10);
            var values = labels.map(l => taxonomyDist[l]);
            
            Plotly.newPlot('taxonomy-pie-plot', [{{
                type: 'pie',
                labels: labels,
                values: values,
                textinfo: 'label+percent',
                hoverinfo: 'label+value'
            }}], {{
                title: 'Gene Taxonomy Distribution (per-gene BLAST)',
                height: 400
            }});
        }}
        
        // Draw visualizations
        if (contigsData.length > 0) {{
            drawCircularGenome();
            drawLinearContigs();
        }}
    </script>
</body>
</html>
"""
    return html


def generate_anomaly_alert_section(results):
    """Generate prominent anomaly alert at top of report"""
    anomaly = results.get('anomaly', {})
    foreign_genes = anomaly.get('gene_taxonomy', {}).get('foreign_genes', [])
    novel_genes = anomaly.get('gene_taxonomy', {}).get('novel_genes', [])
    foreign_clusters = anomaly.get('foreign_gene_clusters', [])
    anomaly_types = anomaly.get('anomaly_types', [])
    severity = anomaly.get('severity', 'none')
    deleted_genes = anomaly.get('alignment', {}).get('deleted_genes', [])
    minority_species = results.get('taxonomy', {}).get('kraken2', {}).get('minority_species', [])
    
    severity_colors = {
        'high': '#c0392b',
        'medium': '#e67e22',
        'low': '#f1c40f'
    }
    severity_color = severity_colors.get(severity, '#e74c3c')
    
    html = f"""
    <div class="section anomaly-alert" style="border-color: {severity_color};">
        <h2 style="color: {severity_color};">⚠️ Genomic Anomalies Detected (Severity: {severity.upper()})</h2>
        <div class="alert alert-danger">
            <strong>This genome shows signs of potential modification or unusual characteristics:</strong>
            <ul style="margin-top: 10px; margin-left: 20px;">
    """
    
    if "foreign_gene_insertion" in anomaly_types:
        html += f"<li><strong>{len(foreign_genes)} foreign gene(s)</strong> detected with different taxonomic origin than host</li>"
    
    if "integrated_foreign_element" in anomaly_types:
        html += f"<li><strong>{len(foreign_clusters)} integrated foreign element(s)</strong> (clusters of 3+ adjacent foreign genes)</li>"
    
    if "novel_synthetic_genes" in anomaly_types:
        html += f"<li><strong>{len(novel_genes)} novel/synthetic gene(s)</strong> with no BLAST hits (potentially engineered)</li>"
    
    if "structural_deletion" in anomaly_types or "gene_deletion" in anomaly_types:
        total_deleted = sum(len(d.get('genes', [])) for d in deleted_genes)
        html += f"<li><strong>Gene deletions</strong> detected - {total_deleted} gene(s) missing relative to reference</li>"
    
    if "gene_rearrangement" in anomaly_types:
        html += "<li><strong>Gene rearrangements</strong> detected - genes out of expected order</li>"
    
    if "gene_inversion" in anomaly_types:
        html += "<li><strong>Gene inversions</strong> detected - genes on opposite strand vs reference</li>"
    
    if "foreign_regulatory_element" in anomaly_types:
        html += "<li><strong>Foreign regulatory elements</strong> detected in intergenic regions</li>"
    
    if "possible_integrated_plasmid" in anomaly_types:
        html += "<li><strong>Possible integrated plasmid/ICE</strong> detected (plasmid backbone genes in chromosome)</li>"
    
    if minority_species:
        html += f"<li><strong>Mixed taxonomy warning:</strong> {len(minority_species)} minority species detected above threshold</li>"
    
    html += """
            </ul>
            <p style="margin-top: 10px;"><em>See detailed analysis sections below for more information.</em></p>
        </div>
    </div>
    """
    return html


def generate_alignment_section(results):
    """Generate reference alignment section if alignment data is available"""
    alignment_regions = results.get('alignment_regions', [])
    anomaly = results.get('anomaly', {})
    ref_info = anomaly.get('reference', {})
    
    if not alignment_regions and ref_info.get('status') != 'found':
        return ""
    
    html = f"""
    <div class="section">
        <h2>Reference Alignment</h2>
        <div class="alert alert-info">
            <strong>Reference:</strong> {ref_info.get('reference', 'N/A')}<br>
            <strong>Mash Distance:</strong> {ref_info.get('distance', 'N/A')}
        </div>
        
        <h3>Alignment Overview</h3>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="value">{len(alignment_regions)}</div>
                <div class="label">Aligned Regions</div>
            </div>
            <div class="stat-card">
                <div class="value">{sum(abs(r.get('ref_end', 0) - r.get('ref_start', 0)) for r in alignment_regions):,}</div>
                <div class="label">Total Aligned (bp)</div>
            </div>
        </div>
        
        <div id="alignment-plot" class="genome-viz-container"></div>
        
        <script>
            // Reference alignment visualization
            if (alignmentRegions.length > 0) {{
                var shapes = [];
                var annotations = [];
                
                // Reference track
                shapes.push({{
                    type: 'rect',
                    x0: 0,
                    x1: Math.max(...alignmentRegions.map(r => r.ref_end)),
                    y0: 1.7,
                    y1: 2.3,
                    fillcolor: '#95a5a6',
                    line: {{ width: 1 }}
                }});
                annotations.push({{ x: -10000, y: 2, text: 'Reference', showarrow: false, xanchor: 'right' }});
                
                // Query track
                shapes.push({{
                    type: 'rect',
                    x0: 0,
                    x1: Math.max(...alignmentRegions.map(r => r.qry_end)),
                    y0: 0.7,
                    y1: 1.3,
                    fillcolor: '#3498db',
                    line: {{ width: 1 }}
                }});
                annotations.push({{ x: -10000, y: 1, text: 'Query', showarrow: false, xanchor: 'right' }});
                
                // Alignment links
                alignmentRegions.forEach(r => {{
                    var color = r.identity > 95 ? '#27ae60' : r.identity > 80 ? '#f39c12' : '#e74c3c';
                    shapes.push({{
                        type: 'path',
                        path: 'M ' + r.ref_start + ' 1.7 L ' + r.ref_end + ' 1.7 L ' + r.qry_end + ' 1.3 L ' + r.qry_start + ' 1.3 Z',
                        fillcolor: color,
                        opacity: 0.5,
                        line: {{ width: 0 }}
                    }});
                }});
                
                Plotly.newPlot('alignment-plot', [], {{
                    shapes: shapes,
                    annotations: annotations,
                    xaxis: {{ title: 'Position (bp)' }},
                    yaxis: {{ visible: false, range: [0, 3] }},
                    title: 'Reference vs Query Alignment',
                    height: 300
                }});
            }}
        </script>
    </div>
    """
    return html


def generate_anomaly_details_section(results):
    """Generate detailed anomaly detection section with all anomaly types"""
    anomaly = results.get('anomaly', {})
    gene_tax = anomaly.get('gene_taxonomy', {})
    foreign_genes = gene_tax.get('foreign_genes', [])
    novel_genes = gene_tax.get('novel_genes', [])
    foreign_clusters = anomaly.get('foreign_gene_clusters', [])
    synteny = results.get('synteny', {})
    intergenic = results.get('intergenic', {})
    integrated = results.get('integrated_elements', {})
    deleted_genes = anomaly.get('alignment', {}).get('deleted_genes', [])
    minority_species = results.get('taxonomy', {}).get('kraken2', {}).get('minority_species', [])
    
    if not anomaly:
        return ""
    
    has_any_anomaly = (foreign_genes or novel_genes or foreign_clusters or 
                       synteny.get('rearrangements') or synteny.get('inversions') or
                       intergenic.get('foreign_elements') or 
                       integrated.get('candidates') or deleted_genes or minority_species)
    
    html = f"""
    <div class="section{' anomaly-alert' if has_any_anomaly else ''}">
        <h2>Anomaly Detection Analysis</h2>
        
        <div class="stats-grid">
            <div class="stat-card{' anomaly' if foreign_genes else ''}">
                <div class="value">{len(foreign_genes)}</div>
                <div class="label">Foreign Genes</div>
            </div>
            <div class="stat-card{' anomaly' if novel_genes else ''}">
                <div class="value">{len(novel_genes)}</div>
                <div class="label">Novel/Unknown Genes</div>
            </div>
            <div class="stat-card{' anomaly' if foreign_clusters else ''}">
                <div class="value">{len(foreign_clusters)}</div>
                <div class="label">Integrated Elements</div>
            </div>
            <div class="stat-card">
                <div class="value">{gene_tax.get('host_genus', 'Unknown')}</div>
                <div class="label">Host Genus</div>
            </div>
        </div>
    """
    
    # Minority Species Alert
    if minority_species:
        html += """
        <h3>⚠️ Minority Species Alert</h3>
        <div class="alert alert-warning">
            <strong>Multiple taxa detected above threshold!</strong> This may indicate contamination, 
            a chimeric/hybrid genome, or horizontal gene transfer.
        </div>
        <table>
            <thead>
                <tr>
                    <th>Genus</th>
                    <th>Percentage</th>
                    <th>Warning</th>
                </tr>
            </thead>
            <tbody>
        """
        for sp in minority_species:
            html += f"""
                <tr style="background: #fff3cd;">
                    <td><strong>{sp.get('name', 'Unknown')}</strong></td>
                    <td>{sp.get('percentage', 0):.2f}%</td>
                    <td><span class="badge badge-warning">Minority Species</span></td>
                </tr>
            """
        html += """
            </tbody>
        </table>
        """
    
    # Foreign Gene Clusters (Integrated Elements)
    if foreign_clusters:
        html += """
        <h3>🧬 Integrated Foreign Elements</h3>
        <div class="alert alert-danger">
            <strong>Clusters of adjacent foreign genes detected!</strong> These may represent 
            integrated plasmids, pathogenicity islands, or engineered gene cassettes.
        </div>
        <table>
            <thead>
                <tr>
                    <th>Element</th>
                    <th>Contig</th>
                    <th>Location</th>
                    <th>Genes</th>
                    <th>Origin</th>
                </tr>
            </thead>
            <tbody>
        """
        for cluster in foreign_clusters:
            html += f"""
                <tr class="foreign-gene">
                    <td><strong>{cluster.get('type', 'Unknown')}</strong></td>
                    <td>{cluster.get('contig', 'N/A')}</td>
                    <td>{cluster.get('start', 0):,} - {cluster.get('end', 0):,}</td>
                    <td>{cluster.get('num_genes', 0)} genes</td>
                    <td><span class="badge badge-danger">{cluster.get('dominant_genus', 'Unknown')}</span></td>
                </tr>
            """
        html += """
            </tbody>
        </table>
        """
    
    # Foreign Genes
    if foreign_genes:
        html += """
        <h3>🔬 Foreign Genes Detected</h3>
        <div class="alert alert-danger">
            The following genes have taxonomic classification different from the host organism.
        </div>
        <table>
            <thead>
                <tr>
                    <th>Gene ID</th>
                    <th>Origin Organism</th>
                    <th>Identity (%)</th>
                    <th>Status</th>
                </tr>
            </thead>
            <tbody>
        """
        for gene in foreign_genes[:20]:
            html += f"""
                <tr class="foreign-gene">
                    <td><strong>{gene.get('gene_id', 'N/A')}</strong></td>
                    <td>{gene.get('organism', 'Unknown')}</td>
                    <td>{gene.get('identity', 0):.1f}%</td>
                    <td><span class="badge badge-danger">Foreign</span></td>
                </tr>
            """
        if len(foreign_genes) > 20:
            html += f'<tr><td colspan="4"><em>... and {len(foreign_genes) - 20} more</em></td></tr>'
        html += """
            </tbody>
        </table>
        """
    
    # Novel/Synthetic Genes (no BLAST hits)
    if novel_genes:
        html += """
        <h3>🧪 Novel/Synthetic Genes</h3>
        <div class="alert alert-danger">
            <strong>Genes with no database matches detected!</strong> These may be synthetic, 
            highly divergent, or from unsequenced organisms.
        </div>
        <table>
            <thead>
                <tr>
                    <th>Gene ID</th>
                    <th>Best Hit</th>
                    <th>Identity</th>
                    <th>Reason</th>
                </tr>
            </thead>
            <tbody>
        """
        for gene in novel_genes[:20]:
            reason = gene.get('reason', 'unknown')
            reason_badge = 'badge-danger' if reason == 'no_blast_hit' else 'badge-warning'
            html += f"""
                <tr style="background: #ffe6e6;">
                    <td><strong>{gene.get('gene_id', 'N/A')}</strong></td>
                    <td>{gene.get('best_hit', 'None') or 'None'}</td>
                    <td>{gene.get('identity', 0):.1f}%</td>
                    <td><span class="badge {reason_badge}">{reason.replace('_', ' ').title()}</span></td>
                </tr>
            """
        if len(novel_genes) > 20:
            html += f'<tr><td colspan="4"><em>... and {len(novel_genes) - 20} more</em></td></tr>'
        html += """
            </tbody>
        </table>
        """
    
    # Deleted Genes
    if deleted_genes:
        html += """
        <h3>❌ Gene Deletions</h3>
        <div class="alert alert-warning">
            <strong>Genes present in reference but absent in this genome.</strong>
        </div>
        <table>
            <thead>
                <tr>
                    <th>Deletion Region</th>
                    <th>Deleted Genes</th>
                </tr>
            </thead>
            <tbody>
        """
        for deletion in deleted_genes[:10]:
            genes_list = ', '.join(g.get('name', 'unknown') for g in deletion.get('genes', [])[:5])
            if len(deletion.get('genes', [])) > 5:
                genes_list += f" ... (+{len(deletion['genes']) - 5} more)"
            html += f"""
                <tr>
                    <td>{deletion.get('deletion_start', 0):,} - {deletion.get('deletion_end', 0):,}</td>
                    <td>{genes_list}</td>
                </tr>
            """
        html += """
            </tbody>
        </table>
        """
    
    # Synteny Analysis
    rearrangements = synteny.get('rearrangements', [])
    inversions = synteny.get('inversions', [])
    if rearrangements or inversions:
        html += f"""
        <h3>🔄 Synteny Analysis</h3>
        <div class="stats-grid">
            <div class="stat-card">
                <div class="value">{synteny.get('total_genes_compared', 0)}</div>
                <div class="label">Genes Compared</div>
            </div>
            <div class="stat-card">
                <div class="value">{synteny.get('genes_in_order', 0)}</div>
                <div class="label">Genes in Order</div>
            </div>
            <div class="stat-card{' anomaly' if rearrangements else ''}">
                <div class="value">{len(rearrangements)}</div>
                <div class="label">Rearrangements</div>
            </div>
            <div class="stat-card{' anomaly' if inversions else ''}">
                <div class="value">{len(inversions)}</div>
                <div class="label">Inversions</div>
            </div>
        </div>
        """
        
        if rearrangements:
            html += """
            <h4>Gene Rearrangements</h4>
            <table>
                <thead><tr><th>Gene</th><th>Type</th></tr></thead>
                <tbody>
            """
            for r in rearrangements[:10]:
                html += f'<tr><td>{r.get("gene", "N/A")}</td><td><span class="badge badge-warning">Out of Order</span></td></tr>'
            html += "</tbody></table>"
        
        if inversions:
            html += """
            <h4>Gene Inversions</h4>
            <table>
                <thead><tr><th>Gene</th><th>Assembly Strand</th><th>Reference Strand</th></tr></thead>
                <tbody>
            """
            for inv in inversions[:10]:
                html += f'<tr><td>{inv.get("gene", "N/A")}</td><td>{inv.get("assembly_strand", "?")}</td><td>{inv.get("ref_strand", "?")}</td></tr>'
            html += "</tbody></table>"
    
    # Intergenic Foreign Elements
    intergenic_foreign = intergenic.get('foreign_elements', [])
    if intergenic_foreign:
        html += """
        <h3>📍 Foreign Intergenic Elements</h3>
        <div class="alert alert-warning">
            <strong>Foreign elements detected in non-coding regions.</strong> 
            These may be regulatory elements (promoters, terminators) from other organisms.
        </div>
        <table>
            <thead>
                <tr>
                    <th>Region</th>
                    <th>Origin</th>
                    <th>Identity</th>
                </tr>
            </thead>
            <tbody>
        """
        for elem in intergenic_foreign[:10]:
            html += f"""
                <tr>
                    <td>{elem.get('region_id', 'N/A')}</td>
                    <td>{elem.get('organism', 'Unknown')}</td>
                    <td>{elem.get('identity', 0):.1f}%</td>
                </tr>
            """
        html += """
            </tbody>
        </table>
        """
    
    # Integrated Mobile Elements (plasmid backbone genes in chromosome)
    integrated_candidates = integrated.get('candidates', []) if isinstance(integrated, dict) else []
    if integrated_candidates:
        html += """
        <h3>🧫 Possible Integrated Plasmids/ICEs</h3>
        <div class="alert alert-warning">
            <strong>Plasmid backbone genes detected in chromosomal contigs.</strong> 
            This suggests integrated mobile elements.
        </div>
        <table>
            <thead>
                <tr>
                    <th>Element</th>
                    <th>Location</th>
                    <th>Backbone Genes</th>
                    <th>Type</th>
                </tr>
            </thead>
            <tbody>
        """
        for elem in integrated_candidates:
            products = ', '.join(elem.get('products', [])[:3])
            html += f"""
                <tr>
                    <td><strong>{elem.get('id', 'N/A')}</strong></td>
                    <td>{elem.get('contig', 'N/A')}: {elem.get('start', 0):,} - {elem.get('end', 0):,}</td>
                    <td>{elem.get('num_backbone_genes', 0)} ({products}...)</td>
                    <td><span class="badge badge-anomaly">{elem.get('type', 'Unknown')}</span></td>
                </tr>
            """
        html += """
            </tbody>
        </table>
        """
    
    if not has_any_anomaly:
        html += """
        <div class="alert alert-info">
            <strong>No significant anomalies detected.</strong> All genes appear consistent with 
            the host organism taxonomy and expected genome structure.
        </div>
        """
    
    html += "</div>"
    return html


def generate_executive_summary(results):
    """Generate executive summary HTML"""
    assembly = results.get('assembly', {})
    taxonomy = results.get('taxonomy', {})
    amr = results.get('amr', {})
    virulence = results.get('virulence', {})
    anomaly = results.get('anomaly', {})
    
    # Get top species
    top_species = "Unknown"
    kraken = taxonomy.get('kraken2', {})
    if kraken.get('species'):
        top_species = kraken['species']
    elif kraken.get('top_hits'):
        top_species = kraken['top_hits'][0].get('name', 'Unknown')
    
    # Anomaly info
    foreign_genes = anomaly.get('gene_taxonomy', {}).get('foreign_genes', [])
    has_anomaly = anomaly.get('anomaly_detected', False)
    
    html = f"""
    <div class="summary-item">
        <span>Organism</span>
        <strong>{top_species}</strong>
    </div>
    <div class="summary-item">
        <span>Genome Size</span>
        <strong>{assembly.get('total_length', 0):,} bp</strong>
    </div>
    <div class="summary-item">
        <span>Contigs</span>
        <strong>{assembly.get('num_contigs', 0)}</strong>
    </div>
    <div class="summary-item">
        <span>N50</span>
        <strong>{assembly.get('n50', 0):,} bp</strong>
    </div>
    <div class="summary-item">
        <span>GC Content</span>
        <strong>{assembly.get('gc_content', 0)}%</strong>
    </div>
    <div class="summary-item">
        <span>AMR Genes</span>
        <strong>{amr.get('summary', {}).get('total_genes', 0)}</strong>
    </div>
    <div class="summary-item">
        <span>Virulence Factors</span>
        <strong>{virulence.get('summary', {}).get('total_factors', 0)}</strong>
    </div>
    <div class="summary-item" style="{'background: rgba(231,76,60,0.3); padding: 8px; border-radius: 4px;' if has_anomaly else ''}">
        <span>{'⚠️ ' if has_anomaly else ''}Foreign Genes</span>
        <strong style="{'color: #e74c3c;' if has_anomaly else ''}">{len(foreign_genes)}{' DETECTED' if has_anomaly else ''}</strong>
    </div>
    """
    return html


def generate_assembly_section(results):
    """Generate assembly section HTML"""
    assembly = results.get('assembly', {})
    
    html = f"""
    <div class="stats-grid">
        <div class="stat-card">
            <div class="value">{assembly.get('num_contigs', 0)}</div>
            <div class="label">Contigs</div>
        </div>
        <div class="stat-card">
            <div class="value">{assembly.get('total_length', 0):,}</div>
            <div class="label">Total Length (bp)</div>
        </div>
        <div class="stat-card">
            <div class="value">{assembly.get('n50', 0):,}</div>
            <div class="label">N50 (bp)</div>
        </div>
        <div class="stat-card">
            <div class="value">{assembly.get('gc_content', 0)}%</div>
            <div class="label">GC Content</div>
        </div>
        <div class="stat-card">
            <div class="value">{assembly.get('largest_contig', 0):,}</div>
            <div class="label">Largest Contig (bp)</div>
        </div>
    </div>
    <p><strong>Assembler:</strong> {assembly.get('assembler', 'N/A')}</p>
    """
    return html


def generate_taxonomy_section(results):
    """Generate taxonomy section HTML with minority species alert"""
    taxonomy = results.get('taxonomy', {})
    kraken = taxonomy.get('kraken2', {})
    
    top_hits = kraken.get('top_hits', [])[:5]
    minority_species = kraken.get('minority_species', [])
    has_minority_alert = kraken.get('minority_species_alert', False)
    
    html = f"""
    <div class="alert alert-info">
        <strong>Organism Type:</strong> {results.get('organism_type', 'Unknown').title()}
    </div>
    """
    
    # Minority species warning
    if has_minority_alert and minority_species:
        html += f"""
        <div class="alert alert-warning">
            <strong>⚠️ Multiple Taxa Detected!</strong> {len(minority_species)} non-dominant genus(es) 
            detected above threshold. This may indicate contamination, chimeric assembly, or 
            extensive horizontal gene transfer.
        </div>
        """
    
    html += """
    <h3>Top Taxonomic Hits</h3>
    <table>
        <thead>
            <tr>
                <th>Organism</th>
                <th>Rank</th>
                <th>Percentage</th>
                <th>Status</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for hit in top_hits:
        rank_name = {'S': 'Species', 'G': 'Genus', 'F': 'Family'}.get(hit.get('rank', ''), hit.get('rank', ''))
        html += f"""
            <tr>
                <td>{hit.get('name', 'Unknown')}</td>
                <td><span class="badge badge-info">{rank_name}</span></td>
                <td>{hit.get('percentage', 0):.1f}%</td>
                <td><span class="badge badge-success">Dominant</span></td>
            </tr>
        """
    
    # Add minority species to table
    for sp in minority_species:
        html += f"""
            <tr style="background: #fff3cd;">
                <td>{sp.get('name', 'Unknown')}</td>
                <td><span class="badge badge-info">Genus</span></td>
                <td>{sp.get('percentage', 0):.2f}%</td>
                <td><span class="badge badge-warning">Minority</span></td>
            </tr>
        """
    
    html += """
        </tbody>
    </table>
    """
    return html


def generate_annotation_section(results):
    """Generate annotation section HTML"""
    annotation = results.get('annotation', {})
    
    html = f"""
    <div class="stats-grid">
        <div class="stat-card">
            <div class="value">{annotation.get('total_genes', 0)}</div>
            <div class="label">Total Genes</div>
        </div>
        <div class="stat-card">
            <div class="value">{annotation.get('cds_count', 0)}</div>
            <div class="label">CDS</div>
        </div>
        <div class="stat-card">
            <div class="value">{annotation.get('trna_count', 0)}</div>
            <div class="label">tRNAs</div>
        </div>
        <div class="stat-card">
            <div class="value">{annotation.get('rrna_count', 0)}</div>
            <div class="label">rRNAs</div>
        </div>
    </div>
    <p><strong>Annotator:</strong> {annotation.get('annotator', 'N/A')}</p>
    """
    return html


def generate_amr_section(results):
    """Generate AMR section HTML"""
    amr = results.get('amr', {})
    genes = amr.get('genes', [])[:10]
    drug_classes = amr.get('drug_classes', [])
    
    if not genes:
        return '<div class="alert alert-info">No antimicrobial resistance genes detected.</div>'
    
    html = f"""
    <div class="alert alert-warning">
        <strong>{len(amr.get('genes', []))} AMR genes detected</strong> across {len(drug_classes)} drug class(es)
    </div>
    
    <h3>Drug Classes</h3>
    <p>{', '.join(drug_classes) if drug_classes else 'None detected'}</p>
    
    <h3>AMR Genes</h3>
    <table>
        <thead>
            <tr>
                <th>Gene</th>
                <th>Drug Class</th>
                <th>Identity</th>
                <th>Coverage</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for gene in genes:
        html += f"""
            <tr>
                <td>{gene.get('gene_symbol', 'N/A')}</td>
                <td>{gene.get('drug_class', 'N/A')}</td>
                <td>{gene.get('identity', 0):.1f}%</td>
                <td>{gene.get('coverage', 0):.1f}%</td>
            </tr>
        """
    
    html += """
        </tbody>
    </table>
    """
    return html


def generate_virulence_section(results):
    """Generate virulence section HTML"""
    virulence = results.get('virulence', {})
    factors = virulence.get('virulence_factors', [])[:10]
    
    if not factors:
        return '<div class="alert alert-info">No virulence factors detected.</div>'
    
    html = f"""
    <div class="alert alert-warning">
        <strong>{virulence.get('summary', {}).get('total_factors', 0)} virulence factors detected</strong>
    </div>
    
    <h3>Virulence Factors</h3>
    <table>
        <thead>
            <tr>
                <th>Factor</th>
                <th>Category</th>
                <th>Identity</th>
            </tr>
        </thead>
        <tbody>
    """
    
    for factor in factors:
        html += f"""
            <tr>
                <td>{factor.get('vfdb_id', 'N/A')}</td>
                <td>{factor.get('category', 'Unknown')}</td>
                <td>{factor.get('identity', 0):.1f}%</td>
            </tr>
        """
    
    html += """
        </tbody>
    </table>
    """
    return html


def generate_mge_section(results):
    """Generate mobile elements section HTML"""
    mge = results.get('mge', {})
    plasmids = mge.get('plasmids', [])
    summary = mge.get('summary', {})
    
    html = f"""
    <div class="stats-grid">
        <div class="stat-card">
            <div class="value">{summary.get('total_plasmids', 0)}</div>
            <div class="label">Plasmids</div>
        </div>
        <div class="stat-card">
            <div class="value">{summary.get('chromosome_contigs', 0)}</div>
            <div class="label">Chromosome Contigs</div>
        </div>
        <div class="stat-card">
            <div class="value">{summary.get('plasmid_contigs', 0)}</div>
            <div class="label">Plasmid Contigs</div>
        </div>
    </div>
    """
    
    if plasmids:
        html += """
        <h3>Identified Plasmids</h3>
        <table>
            <thead>
                <tr>
                    <th>Cluster ID</th>
                    <th>Size (bp)</th>
                    <th>Replicon Type</th>
                    <th>Mobility</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for plasmid in plasmids:
            html += f"""
                <tr>
                    <td>{plasmid.get('cluster_id', 'N/A')}</td>
                    <td>{plasmid.get('size', 0):,}</td>
                    <td>{plasmid.get('rep_type', 'Unknown')}</td>
                    <td>{plasmid.get('predicted_mobility', 'Unknown')}</td>
                </tr>
            """
        
        html += """
            </tbody>
        </table>
        """
    else:
        html += '<div class="alert alert-info">No plasmids detected.</div>'
    
    return html


def generate_prophage_section(results):
    """Generate prophage section HTML"""
    prophage = results.get('prophages', {})
    prophages = prophage.get('prophages', [])
    summary = prophage.get('summary', {})
    
    html = f"""
    <div class="stats-grid">
        <div class="stat-card">
            <div class="value">{summary.get('total_prophages', 0)}</div>
            <div class="label">Prophage Regions</div>
        </div>
        <div class="stat-card">
            <div class="value">{summary.get('total_length', 0):,}</div>
            <div class="label">Total Prophage Length (bp)</div>
        </div>
    </div>
    """
    
    if prophages:
        html += """
        <h3>Prophage Regions</h3>
        <table>
            <thead>
                <tr>
                    <th>ID</th>
                    <th>Contig</th>
                    <th>Start</th>
                    <th>End</th>
                    <th>Length (bp)</th>
                </tr>
            </thead>
            <tbody>
        """
        
        for p in prophages[:10]:
            html += f"""
                <tr>
                    <td>{p.get('id', 'N/A')}</td>
                    <td>{p.get('contig', 'N/A')}</td>
                    <td>{p.get('start', 0):,}</td>
                    <td>{p.get('end', 0):,}</td>
                    <td>{p.get('length', 0):,}</td>
                </tr>
            """
        
        html += """
            </tbody>
        </table>
        """
    else:
        html += '<div class="alert alert-info">No prophage regions detected.</div>'
    
    return html

