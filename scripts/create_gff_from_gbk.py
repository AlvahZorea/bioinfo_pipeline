#!/usr/bin/env python3
"""Create GFF3 from GenBank file"""
from Bio import SeqIO

gff_lines = ['##gff-version 3']
for record in SeqIO.parse('bakta.gbff', 'genbank'):
    for feature in record.features:
        if feature.type in ['CDS', 'gene', 'tRNA', 'rRNA']:
            start = int(feature.location.start) + 1
            end = int(feature.location.end)
            strand = '+' if feature.location.strand == 1 else '-'
            
            attrs = []
            if 'locus_tag' in feature.qualifiers:
                attrs.append(f"locus_tag={feature.qualifiers['locus_tag'][0]}")
                attrs.append(f"Name={feature.qualifiers['locus_tag'][0]}")
            if 'product' in feature.qualifiers:
                # Escape special characters in product
                product = feature.qualifiers['product'][0].replace(';', '%3B').replace('=', '%3D')
                attrs.append(f"product={product}")
            
            gff_lines.append(f"{record.id}\tBakta\t{feature.type}\t{start}\t{end}\t.\t{strand}\t.\t{';'.join(attrs)}")

# Write GFF3 (limit to first 10000 features for performance)
with open('output/citrobacter_sample/04_annotation/bakta/citrobacter_sample.gff3', 'w') as f:
    f.write('\n'.join(gff_lines[:10001]))

print(f'Created GFF3 with {min(len(gff_lines)-1, 10000)} features')

