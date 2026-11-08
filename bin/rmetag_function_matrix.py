#!/usr/bin/env python3
"""Aggregate kallisto estimates using a user-curated gene-to-function mapping."""
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def matrix(path, rows, samples):
    with open(path, 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['feature_id', *samples])
        for feature in sorted(rows):
            writer.writerow([feature, *(f'{rows[feature].get(s, 0):.10g}' for s in samples)])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--mapping', required=True)
    parser.add_argument('--quant', nargs='+', required=True, help='Per-sample kallisto output directories')
    args = parser.parse_args()
    mapping = defaultdict(set)
    with open(args.mapping, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        if not {'gene_id', 'function_id'}.issubset(reader.fieldnames or []):
            parser.error('Mapping requires gene_id and function_id columns')
        for line, row in enumerate(reader, 2):
            gene, function = (row[k].strip() for k in ('gene_id', 'function_id'))
            if not gene or not function:
                parser.error(f'Empty gene_id or function_id at mapping line {line}')
            mapping[gene].add(function)
    if not mapping:
        parser.error('Mapping has no entries')
    genes, functions = defaultdict(dict), defaultdict(lambda: defaultdict(float))
    samples, metrics = [], []
    for directory in args.quant:
        directory = Path(directory)
        sample = directory.name
        if sample in samples:
            parser.error(f'Duplicate sample output: {sample}')
        samples.append(sample)
        total, mapped, found, annotated = 0.0, 0.0, set(), 0
        with open(directory / 'abundance.tsv', newline='') as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            if not {'target_id', 'tpm'}.issubset(reader.fieldnames or []):
                parser.error(f'Invalid kallisto abundance.tsv for {sample}')
            for row in reader:
                gene = row['target_id']
                if gene in found:
                    parser.error(f'Duplicate target {gene} in {sample}')
                found.add(gene)
                value = float(row['tpm'])
                if not math.isfinite(value) or value < 0:
                    parser.error(f'Invalid TPM for {gene} in {sample}')
                genes[gene][sample] = value
                total += value
                if gene in mapping:
                    annotated += 1
                    mapped += value
                    for function in mapping[gene]:
                        functions[function][sample] += value / len(mapping[gene])
        if not found:
            parser.error(f'Empty gene table for {sample}')
        if not annotated:
            parser.error(f'No gene IDs from {sample} occur in the function mapping')
        metrics.append((sample, len(found), annotated, total, mapped))
    samples.sort()
    matrix('gene_tpm.tsv', genes, samples)
    matrix('function_fractional_tpm.tsv', functions, samples)
    with open('annotation_summary.tsv', 'w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['sample', 'genes_in_index', 'genes_with_function', 'total_tpm', 'annotated_tpm'])
        writer.writerows(sorted(metrics))


if __name__ == '__main__':
    main()
