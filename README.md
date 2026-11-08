# RmetaG

RmetaG is a prototype Nextflow workflow for assembly-free characterization of
paired-end shotgun metagenomes. The default is fastp/FastQC, Sylph genome-sketch
profiling and MultiQC. Optional branches add Kraken2/Bracken, host removal,
kallisto gene pseudoalignment, RGI, MetaPhlAn/HUMAnN and DIAMOND.

**Status:** the workflow has passed stub tests, not biological validation on
real samples or production databases. Gene pseudoalignment is experimental for
metagenomic DNA. Outputs must not be used for clinical diagnosis. Database
formats and software versions must be checked before real runs.

## Requirements and samples

Install Nextflow >=25.04 and a supported Java version. Use Docker, Apptainer or
Singularity; software is containerized per process. The test environment used
Nextflow 26.10.5 and Java 21. Production databases and FASTQs are not included.

Input is a CSV with one paired-end sample per row:

```csv
sample,fastq_1,fastq_2
soil_01,/data/soil_01_R1.fastq.gz,/data/soil_01_R2.fastq.gz
soil_02,/data/soil_02_R1.fastq.gz,/data/soil_02_R2.fastq.gz
```

Paths may be absolute or relative to the samplesheet directory. Use gzip FASTQ
and unique IDs beginning with a letter or number and containing only letters,
numbers, dots, underscores and hyphens. Combine lanes before running. The
workflow currently supports paired reads only and uses one reference bundle
for all samples in a run. Run different biome reference configurations
separately. R1/R2 names must be correctly supplied; stub checks do not validate
pairing or FASTQ content.

## Run the fast workflow

```bash
nextflow run main.nf -profile apptainer \
  --input samplesheet.csv \
  --sylph_db /references/soil/soil.syldb \
  --outdir results
```

Use `-profile docker` or `-profile singularity` as appropriate. Without a
container profile, install each tool on `PATH`. Use a Sylph database compatible
with **Sylph 0.9.0**. In particular, `.syl2db` indexes require Sylph >=1.0 and
are incompatible with the included module. Sylph profiles report reference
genome matches; the matrix contains per-reference taxonomic-abundance percent,
not mapped lineage names. Absent rows are represented as zero detections, not
proof of biological absence. Supply and interpret catalogue metadata separately.

## Optional branches

Each optional branch runs only when its parameters are supplied. They consume
the same cleaned reads, after host filtering when enabled.

| Parameters | Branch and output | Notes |
|---|---|---|
| `--host_index DIR --host_index_name NAME` | Hostile 2.0.2 decontamination | Alignment-based. The directory must be a prebuilt Hostile cache; the same host is applied to the entire run. |
| `--kraken_db DIR [--bracken_read_length 150]` | Kraken2 and species-level Bracken | Requires Kraken index and matching Bracken distribution. Set adequate memory. |
| `--gene_index FILE --gene_mapping FILE` | kallisto 0.52.0 gene pseudoalignment and function matrices | Experimental for metagenomic DNA; requires nucleotide CDS references and a curated ID mapping. |
| `--card_db DIR [--rgi_wildcard DIR]` | RGI BWT against CARD | Alignment-based. CARD directory must contain the `card.json` and versioned annotation FASTAs required by RGI. |
| `--metaphlan_db DIR --humann_nucleotide_db DIR --humann_protein_db DIR` | MetaPhlAn 4.1.1 and HUMAnN 3.6.1 | Alignment-based and potentially expensive. Verify software/database compatibility on real data. |
| `--diamond_db FILE.dmnd` | DIAMOND blastx hits, one file per mate | Alignment-based; potentially expensive. Produces raw hits, not functional abundance or pathway tables. |

For example:

```bash
nextflow run main.nf -profile apptainer \
  --input samplesheet.csv --sylph_db /refs/soil.syldb \
  --gene_index /refs/soil-genes.idx --gene_mapping /refs/soil-gene-functions.tsv \
  --outdir results
```

Gene mapping is a tab-separated file with a header and one association per row:

```text
gene_id\tfunction_id
gene_001\tK00001
gene_002\tK00002
```

The example uses `\t` to denote actual tab characters. IDs must match the FASTA
used to build the index. Repeated identical associations are deduplicated. For a
gene mapped to multiple functions its kallisto TPM is divided equally among
those functions; the summed matrix is `function_fractional_tpm.tsv`. The
per-gene TPM matrix and annotation coverage summary are also reported. These
are pseudoalignment-model estimates, not direct gene counts, pathway activities,
or experimentally validated functional abundances. Genes missing from the map
are retained in the gene matrix and excluded from the function matrix.
Redundancy, divergence, multimapping and fragment modeling need benchmarking.

HUMAnN receives a single gzip FASTQ made by concatenating both clean mates as
separate gzip members; MetaPhlAn profiles and HUMAnN inputs are joined by sample
ID. Confirm compatible MetaPhlAn/ChocoPhlAn releases and verify handling of your
read identifiers on real data. Its per-sample gene-family and pathway files are
published in native units; cohort HUMAnN normalization is not yet implemented.
The DIAMOND branch publishes raw alignments only; it does not aggregate KOs.
BioSIFTR-style imputed functions are not included.

## Database helper

`scripts/databases.py` requires Python 3.9+ and has no third-party Python
dependencies. It never automatically downloads large databases during a pipeline
run. To list the third-party MGnify-derived Sylph indexes, built with Sylph
0.9.0 and published in [Zenodo record 21842023](https://doi.org/10.5281/zenodo.21842023):

```bash
scripts/databases.py list-mgnify
scripts/databases.py fetch-mgnify --biome soil --outdir /refs/sylph
scripts/databases.py fetch-mgnify --biome human-gut --outdir /refs/sylph
```

The soil file is approximately 3 GB, the human-gut file approximately 454 MB;
check available disk space and the catalogue's source/release for your study.
The helper retrieves the record's published MD5, verifies the file, and writes
a JSON manifest with its SHA-256, URL and timestamp. It will not overwrite an
existing destination. These biome indexes are third-party derivatives of MGnify
catalogues, not a comprehensive all-domain reference or automatically matched
taxonomy metadata.

For other references, download an HTTPS artifact with a publisher-supplied
checksum; archives are **not** extracted automatically:

```bash
scripts/databases.py fetch-url --url https://example.org/reference.tar.gz \
  --checksum sha256:REPLACE_WITH_PUBLISHED_HEX --output /refs/reference.tar.gz
```

Use an actual URL and published checksum. For Kraken2/Bracken, choose a pinned
[reference package](https://benlangmead.github.io/aws-indexes/k2), inspect its
contents and Bracken read-length distributions, download it with the published
checksum, and extract it separately into a dedicated directory. CARD and
WildCARD resources must likewise be prepared according to their own download
and licensing instructions; `fetch-url` does not construct their RGI annotation
files. User-selected FASTA and function mappings must be curated and compatible.

Build a gene index with the same kallisto version used by the workflow:

```bash
scripts/databases.py index-genes --fasta /refs/coding-sequences.fna.gz \
  --output /refs/genes.idx
```

This invokes the installed `kallisto`, records its reported version and the
FASTA SHA-256, and does not generate the mapping file. Use nucleotide CDS, not
protein sequences. Index creation and the chosen references need validation.

The helper can also call installed upstream database tools. Match their versions
to the pipeline containers and choose compatible, explicitly named releases:

```bash
scripts/databases.py fetch-host --name human-t2t-hla --outdir /refs/hostile
scripts/databases.py fetch-metaphlan --index EXACT_INDEX_NAME --outdir /refs/metaphlan
scripts/databases.py fetch-humann --database chocophlan --build full --outdir /refs/chocophlan
scripts/databases.py fetch-humann --database uniref --build uniref90_diamond --outdir /refs/uniref
```

These wrappers require the respective executables on `PATH`; they delegate
downloads to the upstream tools and do not themselves verify upstream tool
downloads or establish cross-database compatibility. HUMAnN is invoked with
`--update-config no`. Reference versions and licenses must be recorded by the
user. The wrappers have not been tested with full reference downloads.

## Results and testing

Outputs include `qc/`, `taxonomy/`, `function/`, `resistance/`, `matrices/`,
and `multiqc/` according to enabled branches. Cleaned FASTQs remain in the
Nextflow work directory rather than being published; retain it for `-resume`.
Configure memory, CPUs and storage for the actual indexes and 30 Gb samples.

Check basic and all optional branch wiring with tiny fake input and reference
files:

```bash
tests/run-stub.sh
python3 -m unittest discover -s tests -p 'test_*.py'
```

The stub runner requires Nextflow and Python, uses temporary work/output
folders, and accepts `NEXTFLOW_BIN` to select a launcher. Stub tests do not run
the real tools and establish neither database compatibility nor accuracy.
The real helper was tested against Zenodo with a small human-vaginal index;
large soil/reference downloads and upstream downloader wrappers have not been
tested. The next milestone is a real-data benchmark, including excluded-reference
controls and comparison to HUMAnN/DIAMOND on representative subsets.

Selected nf-core modules and modifications are documented in
[modules/README.md](modules/README.md).

The repository's `.gitignore` excludes local Nextflow runs, raw FASTQs and
reference indexes while retaining the small fixtures under `tests/`. Keep
production sample data and downloaded reference bundles outside the repository.
