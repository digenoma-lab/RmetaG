process CONCAT_READS {
    tag "$meta.id"
    container 'ubuntu:24.04'
    input:
    tuple val(meta), path(reads)
    output:
    tuple val(meta), path("${meta.id}.combined.fastq.gz"), emit: reads
    script:
    """
    cat '${reads[0]}' '${reads[1]}' > '${meta.id}.combined.fastq.gz'
    """
}
