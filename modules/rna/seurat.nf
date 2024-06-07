/*
Normalizing and clustering data with seurat
*/

params.options = [:]

process SEURAT {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/seurat/", pattern: '*.rds', mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_xlarge'
    beforeScript 'ulimit -u $(ulimit -Hu)'

    input:
    tuple val(sampleId), path(matrix), path(barcodes), path(features)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}_seurat.rds"), emit: rds
    tuple val(sampleId), path("${sampleId}_umap.csv"), emit: umap_csv
    tuple val(sampleId), path("${sampleId}_top_features.csv"), emit: top_features_csv
    tuple val(sampleId), path("${sampleId}_seurat_metadata.csv.gz"), emit: metadata
    path("${sampleId}_SEURAT_messages.txt"), optional: true, emit: messages

    script:
    """
    rnaSeurat.R ${sampleId} ./ ${params.options.mixed}
    """
}
