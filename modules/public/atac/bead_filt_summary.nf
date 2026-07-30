/*
Bead filtration summary table
*/

process BEAD_FILT_SUMMARY {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(stats), path(quant), path(merge_params)
    path counts, stageAs: 'fastqTIreadcounts.csv'
    val images_pulled

    output:
    path("${sampleId}.beadFiltSummary.csv"), emit:summary

    script:
    """
    publicAtacBeadFiltSummary.py -s ${sampleId} -i ./
    """
}
