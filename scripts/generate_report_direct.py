#!/usr/bin/env python3
"""
Generate the HTML report directly (bypassing Snakemake)
for testing visualization with existing outputs.
"""

import json
import os
import sys
from pathlib import Path
from Bio import SeqIO

# Add workflow path for imports
sys.path.insert(0, str(Path(__file__).parent.parent))

SAMPLE_ID = "citrobacter_sample"
OUTPUT_DIR = Path(f"output/{SAMPLE_ID}")

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


def load_all_results():
    """Load all results from JSON files"""
    results = {
        "sample_id": SAMPLE_ID,
        "report_generated": "2026-01-01T13:00:00",
        "pipeline_version": "1.0.0",
        "organism_type": "bacteria"
    }
    
    # Load each JSON file if it exists
    json_files = {
        "assembly": OUTPUT_DIR / "02_assembly" / "assembly_stats.json",
        "taxonomy": OUTPUT_DIR / "03_identification" / "taxonomy.json",
        "annotation": OUTPUT_DIR / "04_annotation" / "annotation_stats.json",
        "amr": OUTPUT_DIR / "05_specialized" / "amr" / "amr_summary.json",
        "virulence": OUTPUT_DIR / "05_specialized" / "virulence" / "virulence_summary.json",
        "mge": OUTPUT_DIR / "05_specialized" / "mge" / "mge_summary.json",
        "prophages": OUTPUT_DIR / "05_specialized" / "prophages" / "prophage_summary.json",
        "anomaly": OUTPUT_DIR / "06_anomaly" / "anomaly_summary.json",
    }
    
    for key, filepath in json_files.items():
        if filepath.exists():
            with open(filepath) as f:
                results[key] = json.load(f)
            print(f"  Loaded: {key}")
        else:
            print(f"  Missing: {key}")
    
    # Load contigs
    contigs_file = OUTPUT_DIR / "02_assembly" / "contigs.fasta"
    contigs_data = []
    for record in SeqIO.parse(contigs_file, "fasta"):
        seq = str(record.seq).upper()
        gc = (seq.count('G') + seq.count('C')) / len(seq) * 100 if seq else 0
        contigs_data.append({
            "id": record.id,
            "length": len(record),
            "gc": round(gc, 2)
        })
    results["contigs"] = contigs_data
    
    # Load features from GFF
    gff_file = OUTPUT_DIR / "04_annotation" / "bakta" / f"{SAMPLE_ID}.gff3"
    if gff_file.exists():
        results["features"] = parse_gff_features(gff_file)[:500]  # Limit for performance
        print(f"  Loaded: {len(results['features'])} features")
    else:
        results["features"] = []
    
    results["alignment_regions"] = []  # No reference alignment in this test
    
    return results


def generate_html(results):
    """Generate the HTML report"""
    contigs_data = results.get("contigs", [])
    features = results.get("features", [])
    alignment_regions = results.get("alignment_regions", [])
    anomaly = results.get("anomaly", {})
    
    has_anomaly = anomaly.get("anomaly_detected", False)
    foreign_genes = anomaly.get("gene_taxonomy", {}).get("foreign_genes", [])
    
    contig_lengths = [c["length"] for c in contigs_data]
    
    title = f"Genome Analysis Report: {SAMPLE_ID}"
    
    # Generate sections
    assembly = results.get('assembly', {})
    taxonomy = results.get('taxonomy', {})
    annotation = results.get('annotation', {})
    
    # Get top species
    top_species = "Unknown"
    kraken = taxonomy.get('kraken2', {})
    if kraken.get('species'):
        top_species = kraken['species']
    elif kraken.get('top_hits'):
        top_species = kraken['top_hits'][0].get('name', 'Unknown')
    
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
                <strong>Sample ID:</strong> {SAMPLE_ID} | 
                <strong>Organism Type:</strong> {results.get('organism_type', 'Unknown')} | 
                <strong>Generated:</strong> {results.get('report_generated', 'N/A')[:19]}
                {"| <span style='color: #ffcccc;'>⚠️ ANOMALIES DETECTED</span>" if has_anomaly else ""}
            </div>
        </header>
        
        {"" if not has_anomaly else f'''
        <div class="section anomaly-alert">
            <h2>⚠️ Genomic Anomalies Detected</h2>
            <div class="alert alert-danger">
                <strong>This genome shows signs of potential modification:</strong>
                <ul style="margin-top: 10px; margin-left: 20px;">
                    <li><strong>{len(foreign_genes)} foreign gene(s)</strong> detected with different taxonomic origin than host</li>
                </ul>
                <p style="margin-top: 10px;"><em>See detailed analysis sections below.</em></p>
            </div>
        </div>
        '''}
        
        <!-- Executive Summary -->
        <div class="section">
            <h2>Executive Summary</h2>
            <div class="summary-box{' alert-box' if has_anomaly else ''}">
                <h3>Key Findings</h3>
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
                <div class="summary-item" style="{'background: rgba(231,76,60,0.3); padding: 8px; border-radius: 4px;' if has_anomaly else ''}">
                    <span>{'⚠️ ' if has_anomaly else ''}Foreign Genes</span>
                    <strong style="{'color: #e74c3c;' if has_anomaly else ''}">{len(foreign_genes)}{' DETECTED' if has_anomaly else ''}</strong>
                </div>
            </div>
        </div>
        
        <!-- Genome Visualization -->
        <div class="section">
            <h2>Genome Visualization</h2>
            <p>Interactive visualization of genome structure and anomalies.</p>
            
            <h3>Circular Genome Map</h3>
            <div id="circular-genome-plot" class="genome-viz-container"></div>
            
            <div class="legend">
                <div class="legend-item"><div class="legend-color" style="background: #3498db;"></div> CDS</div>
                <div class="legend-item"><div class="legend-color" style="background: #27ae60;"></div> tRNA</div>
                <div class="legend-item"><div class="legend-color" style="background: #f39c12;"></div> rRNA</div>
                <div class="legend-item"><div class="legend-color" style="background: #e74c3c;"></div> Foreign Gene</div>
            </div>
            
            <h3>Linear Contig View</h3>
            <div id="linear-contig-plot" class="genome-viz-container"></div>
        </div>
        
        <!-- Anomaly Detection -->
        <div class="section{' anomaly-alert' if foreign_genes else ''}">
            <h2>Anomaly Detection Analysis</h2>
            
            <div class="stats-grid">
                <div class="stat-card{' anomaly' if foreign_genes else ''}">
                    <div class="value">{len(foreign_genes)}</div>
                    <div class="label">Foreign Genes</div>
                </div>
                <div class="stat-card">
                    <div class="value">{anomaly.get('gene_taxonomy', {}).get('host_genus', 'Unknown')}</div>
                    <div class="label">Host Genus</div>
                </div>
            </div>
            
            {"" if not foreign_genes else '''
            <h3>Foreign Genes Detected</h3>
            <div class="alert alert-danger">
                The following genes have taxonomic classification different from the host organism.
            </div>
            <table>
                <thead>
                    <tr><th>Gene ID</th><th>Origin Organism</th><th>Identity (%)</th><th>Status</th></tr>
                </thead>
                <tbody>
            ''' + ''.join(f'''
                    <tr class="foreign-gene">
                        <td><strong>{gene.get('gene_id', 'N/A')}</strong></td>
                        <td>{gene.get('organism', 'Unknown')}</td>
                        <td>{gene.get('identity', 0):.1f}%</td>
                        <td><span class="badge badge-danger">Foreign</span></td>
                    </tr>
            ''' for gene in foreign_genes) + '''
                </tbody>
            </table>
            '''}
        </div>
        
        <!-- Assembly -->
        <div class="section">
            <h2>Assembly Statistics</h2>
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
            </div>
            <div id="contig-length-plot" class="plot-container"></div>
        </div>
        
        <!-- Taxonomy -->
        <div class="section">
            <h2>Taxonomic Identification</h2>
            <div class="alert alert-info">
                <strong>Organism Type:</strong> {results.get('organism_type', 'Unknown').title()}
            </div>
            <h3>Taxonomy Distribution (per-gene)</h3>
            <div id="taxonomy-pie-plot" class="plot-container"></div>
        </div>
        
        <!-- Annotation -->
        <div class="section">
            <h2>Genome Annotation</h2>
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
        </div>
        
        <footer>
            Genome Analysis Pipeline v1.0.0 | Report generated: {results.get('report_generated', 'N/A')[:19]}
        </footer>
    </div>
    
    <script>
        var contigsData = {json.dumps(contigs_data)};
        var features = {json.dumps(features[:200])};
        var foreignGenes = {json.dumps(foreign_genes)};
        var taxonomyDist = {json.dumps(anomaly.get('gene_taxonomy', {}).get('taxonomy_distribution', {}))};
        var contigLengths = {json.dumps(contig_lengths)};
        
        // Circular Genome Map
        function drawCircularGenome() {{
            var totalLength = contigsData.reduce((sum, c) => sum + c.length, 0);
            var traces = [];
            
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
            
            // Foreign genes as markers
            if (foreignGenes.length > 0) {{
                var fTheta = [];
                var fR = [];
                var fText = [];
                
                foreignGenes.forEach((fg, i) => {{
                    // Place markers evenly around the genome
                    fTheta.push((i / foreignGenes.length) * 360);
                    fR.push(0.75);
                    fText.push(fg.gene_id + '<br>Origin: ' + fg.organism + '<br>Identity: ' + fg.identity.toFixed(1) + '%');
                }});
                
                traces.push({{
                    type: 'scatterpolar',
                    mode: 'markers',
                    r: fR,
                    theta: fTheta,
                    marker: {{ size: 12, color: '#e74c3c', symbol: 'diamond' }},
                    text: fText,
                    hoverinfo: 'text',
                    name: 'Foreign Genes'
                }});
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
            var shapes = [];
            var annotations = [];
            
            contigsData.forEach((contig, i) => {{
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
                    x: -50000,
                    y: i,
                    text: contig.id,
                    showarrow: false,
                    xanchor: 'right'
                }});
            }});
            
            var layout = {{
                shapes: shapes,
                annotations: annotations,
                xaxis: {{ title: 'Position (bp)', rangemode: 'tozero' }},
                yaxis: {{ 
                    tickvals: contigsData.map((c, i) => i),
                    ticktext: contigsData.map(c => ''),
                    range: [-1, contigsData.length]
                }},
                title: 'Linear Contig View',
                height: Math.max(300, contigsData.length * 50 + 100),
                showlegend: false
            }};
            
            Plotly.newPlot('linear-contig-plot', [], layout);
        }}
        
        // Contig length histogram
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
        
        if (contigsData.length > 0) {{
            drawCircularGenome();
            drawLinearContigs();
        }}
    </script>
</body>
</html>
"""
    return html


def main():
    print("=" * 60)
    print("Generating report directly (bypassing Snakemake)")
    print("=" * 60)
    
    print("\nLoading results...")
    results = load_all_results()
    
    print("\nGenerating HTML...")
    html = generate_html(results)
    
    output_file = OUTPUT_DIR / "08_report" / "genome_report.html"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html)
    
    print(f"\nReport saved to: {output_file}")
    print("=" * 60)


if __name__ == "__main__":
    main()

