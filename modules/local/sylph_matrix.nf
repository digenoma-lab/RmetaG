process SYLPH_MATRIX {
    tag 'cohort'
    container 'python:3.11-slim'
    input:
    path profiles, stageAs: '?/*'
    output:
    path 'sylph_taxonomic_abundance_percent.tsv', emit: matrix
    script:
    """
    rmetag_sylph_matrix.py ${profiles.collect { profile -> "'${profile}'" }.join(' ')}
    """
}
