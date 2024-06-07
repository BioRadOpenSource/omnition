/*
Label cells as either true or false
*/

params.options = [:]

process CALCULATE_BEADS_PER_DROP {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
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
    atacCalculateBeadsPerDrop.R -o ./  -p ${sampleId} ${qc}
    """
}
