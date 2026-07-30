/*
Create the MAS compatible count matrix
*/

process MAS_MATRIX_CYTO {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/masMatrix/",
        mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_large'

    input:
    tuple val(sampleId), path(matrix), path(barcodes), path(antibodies)
    val images_pulled

    output:
    tuple val(sampleId), path("Matrix.mtx.gz"), path("Barcodes.tsv.gz"),
        path("Features.tsv.gz"), emit: mas_matrix

    script:
    """
    publicCytometryMasOutput.py
    """
}
