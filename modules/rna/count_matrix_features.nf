/*
Counting expression matrix features
*/

params.options = [:]

process COUNT_MATRIX_FEATURES {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(matrix), path(barcodes), path(features)
    val images_pulled

    output:
    tuple val(sampleId), path('*_matrix_features.csv'), emit: count

    script:
    """
    rnaCountMatrixFeatures.R ${sampleId} ${matrix}
    """
}
