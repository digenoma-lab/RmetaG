nextflow.enable.dsl = 2

include { FASTP } from './modules/nf-core/fastp/main'
include { FASTQC as FASTQC_RAW } from './modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CLEAN } from './modules/nf-core/fastqc/main'
include { HOSTILE_CLEAN } from './modules/nf-core/hostile/clean/main'
include { SYLPH_PROFILE } from './modules/nf-core/sylph/profile/main'
include { KRAKEN2_KRAKEN2 } from './modules/nf-core/kraken2/kraken2/main'
include { BRACKEN_BRACKEN } from './modules/nf-core/bracken/bracken/main'
include { MULTIQC } from './modules/nf-core/multiqc/main'
include { KALLISTO_QUANT } from './modules/nf-core/kallisto/quant/main'
include { RGI_BWT } from './modules/nf-core/rgi/bwt/main'
include { FUNCTION_MATRIX } from './modules/local/function_matrix'
include { CONCAT_READS } from './modules/local/concat_reads'
include { METAPHLAN_METAPHLAN } from './modules/nf-core/metaphlan/metaphlan/main'
include { HUMANN3_HUMANN } from './modules/nf-core/humann3/humann/main'
include { DIAMOND_BLASTX } from './modules/nf-core/diamond/blastx/main'
include { SYLPH_MATRIX } from './modules/local/sylph_matrix'

workflow {
    if (params.help) {
        log.info 'RmetaG: --input samples.csv --sylph_db PATH [--gene_index FILE --gene_mapping FILE] [--host_index DIR --host_index_name NAME] [--kraken_db DIR] [--card_db DIR] [--metaphlan_db DIR --humann_nucleotide_db DIR --humann_protein_db DIR] [--diamond_db FILE] --outdir results. See README.md.'
        exit 0
    }
    if (!params.input || !params.sylph_db) error 'Both --input and --sylph_db are required. See README.md.'
    if (!!params.host_index != !!params.host_index_name) error '--host_index and --host_index_name must be provided together.'
    if (!!params.gene_index != !!params.gene_mapping) error '--gene_index and --gene_mapping must be provided together.'
    if (params.rgi_wildcard && !params.card_db) error '--rgi_wildcard requires --card_db.'
    def humannSettings = [params.metaphlan_db, params.humann_nucleotide_db, params.humann_protein_db]
    if (humannSettings.any { setting -> setting } && !humannSettings.every { setting -> setting }) error 'HUMAnN requires --metaphlan_db, --humann_nucleotide_db, and --humann_protein_db.'
    if (params.kraken_db && !(params.bracken_read_length.toString() ==~ /[1-9][0-9]*/)) error '--bracken_read_length must be a positive integer.'
    def sheet = file(params.input, checkIfExists: true)
    def sylphDb = file(params.sylph_db, checkIfExists: true)
    def hostDb = params.host_index ? file(params.host_index, checkIfExists: true) : null
    def krakenDb = params.kraken_db ? file(params.kraken_db, checkIfExists: true) : null
    def geneIndex = params.gene_index ? file(params.gene_index, checkIfExists: true) : null
    def geneMapping = params.gene_mapping ? file(params.gene_mapping, checkIfExists: true) : null
    def cardDb = params.card_db ? file(params.card_db, checkIfExists: true) : null
    def wildcardDb = params.rgi_wildcard ? file(params.rgi_wildcard, checkIfExists: true) : []
    def metaphlanDb = params.metaphlan_db ? file(params.metaphlan_db, checkIfExists: true) : null
    def humannNucleotideDb = params.humann_nucleotide_db ? file(params.humann_nucleotide_db, checkIfExists: true) : null
    def humannProteinDb = params.humann_protein_db ? file(params.humann_protein_db, checkIfExists: true) : null
    def diamondDb = params.diamond_db ? file(params.diamond_db, checkIfExists: true) : null

    ch_reads = channel.fromPath(sheet)
        .splitCsv(header: true, strip: true)
        .toList()
        .flatMap { rows ->
            if (!rows) error 'Samplesheet is empty.'
            def seen = [] as Set
            rows.collect { row ->
                def id = row.sample?.trim()
                if (!id || !(id ==~ /[A-Za-z0-9][A-Za-z0-9._-]*/)) error "Invalid sample ID: ${id}; use letters, digits, dot, underscore or hyphen."
                if (!seen.add(id)) error "Duplicate sample ID: ${id}. Combine lanes before running."
                if (!row.fastq_1?.trim() || !row.fastq_2?.trim()) error "Sample ${id} requires both fastq_1 and fastq_2."
                def read1 = file(row.fastq_1.startsWith('/') ? row.fastq_1 : sheet.parent.resolve(row.fastq_1), checkIfExists: true)
                def read2 = file(row.fastq_2.startsWith('/') ? row.fastq_2 : sheet.parent.resolve(row.fastq_2), checkIfExists: true)
                if (read1 == read2) error "Sample ${id} has the same file for both mates."
                tuple([id: id, single_end: false], [read1, read2])
            }
        }

    FASTQC_RAW(ch_reads)
    FASTP(ch_reads.map { meta, reads -> tuple(meta, reads, []) }, false, false, false)

    ch_host_reports = channel.empty()
    if (hostDb) {
        HOSTILE_CLEAN(FASTP.out.reads, channel.value(tuple(params.host_index_name, hostDb)))
        ch_clean = HOSTILE_CLEAN.out.fastq.map { meta, reads -> tuple(meta, [reads].flatten().sort { read -> read.name }) }
        ch_host_reports = HOSTILE_CLEAN.out.json.map { _meta, report -> report }
    } else {
        ch_clean = FASTP.out.reads
    }

    FASTQC_CLEAN(ch_clean)
    SYLPH_PROFILE(ch_clean, channel.value(sylphDb))
    SYLPH_MATRIX(SYLPH_PROFILE.out.profile_out.map { _meta, profile -> profile }.collect())
    if (krakenDb) {
        KRAKEN2_KRAKEN2(ch_clean, channel.value(krakenDb), false, false)
        BRACKEN_BRACKEN(KRAKEN2_KRAKEN2.out.report, channel.value(krakenDb))
    }

    if (geneIndex) {
        KALLISTO_QUANT(ch_clean, channel.value(tuple([id: 'genes'], geneIndex, [], [])), 0, 0)
        FUNCTION_MATRIX(KALLISTO_QUANT.out.results.map { _meta, result -> result }.collect(), channel.value(geneMapping))
    }
    if (cardDb) {
        RGI_BWT(ch_clean, channel.value(cardDb), channel.value(wildcardDb))
    }
    if (metaphlanDb) {
        METAPHLAN_METAPHLAN(ch_clean, channel.value(metaphlanDb), false)
        CONCAT_READS(ch_clean)
        ch_humann = CONCAT_READS.out.reads.join(METAPHLAN_METAPHLAN.out.profile).multiMap { meta, reads, profile ->
            reads: tuple(meta, reads)
            profile: tuple(meta, profile)
        }
        HUMANN3_HUMANN(ch_humann.reads, ch_humann.profile, channel.value(humannNucleotideDb), channel.value(humannProteinDb), channel.value([]))
    }
    if (diamondDb) {
        ch_mates = ch_clean.flatMap { meta, reads ->
            [tuple(meta + [id: "${meta.id}_R1"], reads[0]), tuple(meta + [id: "${meta.id}_R2"], reads[1])]
        }
        DIAMOND_BLASTX(ch_mates, channel.value(tuple([id: 'proteins'], diamondDb)), 'txt', '')
    }

    ch_reports = FASTQC_RAW.out.zip.map { _meta, files -> files }
        .mix(FASTQC_CLEAN.out.zip.map { _meta, files -> files })
        .mix(FASTP.out.json.map { _meta, report -> report })
        .mix(ch_host_reports)
        .collect(flat: true)
        .map { files -> tuple([id: 'multiqc'], files, [], [], [], []) }
    MULTIQC(ch_reports)
}
