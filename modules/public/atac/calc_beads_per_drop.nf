/*
Label cells as either true or false
*/

process CALCULATE_BEADS_PER_DROP {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(qc)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.cell_data.csv"), emit: barcode
    tuple val(sampleId), path("${sampleId}.cell_data.csv"), emit: cells
    val(sampleId), emit: sampleid

    script:
    """
    publicAtacCalculateBeadsPerDrop.R -o ./  -p ${sampleId} ${qc}
    """
}
