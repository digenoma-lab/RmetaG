#!/usr/bin/env python3
"""Combine per-sample Sylph estimates by reference genome."""
import argparse
import csv
import math
from pathlib import Path
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('profiles', nargs='+')
args = parser.parse_args()
samples = []
values = defaultdict(dict)
for filename in args.profiles:
    sample = Path(filename).stem
    if sample in samples:
        parser.error(f'duplicate sample: {sample}')
    samples.append(sample)
    with open(filename, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'Genome_file', 'Taxonomic_abundance'}.issubset(reader.fieldnames or []):
            parser.error(f'invalid Sylph header in {filename}')
        for row in reader:
            ref = row['Genome_file']
            abundance = float(row['Taxonomic_abundance'])
            if not ref or not math.isfinite(abundance) or abundance < 0:
                parser.error(f'invalid abundance or genome in {filename}')
            if sample in values[ref]:
                parser.error(f'duplicate genome {ref} in {filename}')
            values[ref][sample] = abundance
samples.sort()
with open('sylph_taxonomic_abundance_percent.tsv', 'w', newline='') as handle:
    writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
    writer.writerow(['genome_reference', *samples])
    for ref in sorted(values):
        writer.writerow([ref, *(f'{values[ref].get(sample, 0):.10g}' for sample in samples)])
