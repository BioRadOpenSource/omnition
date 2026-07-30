/*
Create the MAS compatible count matrix
*/

process MAS_MATRIX_CITE {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/masMatrix/",
        mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_large'

    input:
    tuple val(sampleId),
        path(matrix, name: 'cyto_matrix/*'),
        path(barcodes, name: 'cyto_matrix/*'), // groovylint-disable-line
        path(antibodies, name: 'cyto_matrix/*'), // groovylint-disable-line
        path(rna_matrix, name: 'rna_matrix/*'),
        path(rna_barcodes, name: 'rna_matrix/*'), // groovylint-disable-line
        path(rna_genes, name: 'rna_matrix/*') // groovylint-disable-line
    val images_pulled

    output:
    tuple val(sampleId), path("Matrix.mtx.gz"), path("Barcodes.tsv.gz"),
        path("Features.tsv.gz"), emit: mas_matrix

    script:
    """
    publicCytometryMasOutput.py --input-dir ./cyto_matrix --rna-dir ./rna_matrix
    """
}
