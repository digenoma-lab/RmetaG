#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd "$(dirname "$0")/.." && pwd)"
nextflow_bin="${NEXTFLOW_BIN:-nextflow}"
test_dir="$(mktemp -d "${TMPDIR:-/tmp}/rmetag-stub.XXXXXXXX")"
trap 'rm -rf "$test_dir"' EXIT
cd "$repo_dir"
"$nextflow_bin" run main.nf -profile stub -stub-run -w "$test_dir/work-basic" --outdir "$test_dir/basic"
"$nextflow_bin" run main.nf -profile stub_all -stub-run -w "$test_dir/work-all" --outdir "$test_dir/all"
for sample in sample_a sample_b; do
    test -f "$test_dir/basic/taxonomy/sylph/${sample}.tsv"
    test -f "$test_dir/all/taxonomy/sylph/${sample}.tsv"
    test -f "$test_dir/all/taxonomy/kraken2/${sample}.kraken2.report.txt"
    test -f "$test_dir/all/taxonomy/bracken/${sample}.tsv"
    test -f "$test_dir/all/qc/host/${sample}.json"
    test -f "$test_dir/basic/qc/fastp/${sample}.fastp.json"
    test -f "$test_dir/all/function/kallisto/${sample}/abundance.tsv"
    test -f "$test_dir/all/resistance/rgi/${sample}.json"
    test -f "$test_dir/all/function/humann/${sample}_genefamilies.tsv.gz"
    test -f "$test_dir/all/taxonomy/metaphlan/${sample}_profile.txt"
    test -f "$test_dir/all/function/diamond/${sample}_R1.txt"
    test -f "$test_dir/all/function/diamond/${sample}_R2.txt"
done
test -f "$test_dir/basic/multiqc/multiqc_report.html"
test -f "$test_dir/all/multiqc/multiqc_report.html"
test -f "$test_dir/all/matrices/function/gene_tpm.tsv"
test -f "$test_dir/all/matrices/function/function_fractional_tpm.tsv"
test -f "$test_dir/basic/matrices/taxonomy/sylph_taxonomic_abundance_percent.tsv"
python3 - "$test_dir/all/matrices/function/function_fractional_tpm.tsv" <<'PY'
import csv, sys
with open(sys.argv[1]) as f:
    rows = list(csv.reader(f, delimiter='\t'))
assert rows == [['feature_id', 'sample_a', 'sample_b'], ['K00001', '1000000', '1000000']], rows
PY
printf 'sample,fastq_1,fastq_2\ninvalid,missing.fastq.gz,\n' > "$test_dir/invalid.csv"
if "$nextflow_bin" run main.nf -profile stub -stub-run -w "$test_dir/work-invalid" --input "$test_dir/invalid.csv" --outdir "$test_dir/invalid" > "$test_dir/invalid.log" 2>&1; then
    printf 'Invalid samplesheet unexpectedly succeeded\n' >&2
    exit 1
fi
if ! grep -q 'requires both fastq_1 and fastq_2' "$test_dir/invalid.log"; then
    cat "$test_dir/invalid.log" >&2
    exit 1
fi
printf 'Stub workflow checks passed: basic, all optional branches, and invalid input.\n'
