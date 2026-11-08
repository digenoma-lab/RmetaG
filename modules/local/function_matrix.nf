process FUNCTION_MATRIX {
    tag 'cohort'
    container 'python:3.11-slim'

    input:
    path quant_dirs, stageAs: '?/*'
    path mapping

    output:
    path 'gene_tpm.tsv', emit: genes
    path 'function_fractional_tpm.tsv', emit: functions
    path 'annotation_summary.tsv', emit: summary

    script:
    """
    rmetag_function_matrix.py --mapping '${mapping}' --quant ${quant_dirs.collect { dir -> "'${dir}'" }.join(' ')}
    """
}
