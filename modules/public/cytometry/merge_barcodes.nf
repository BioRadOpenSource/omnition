/*
Parsing and correcting cell barcodes
*/

process MERGE_BARCODES {
    tag "${sampleId}"
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/countMatrix/", mode: 'copy',
        overwrite: true, pattern: '*unfiltered.*'
    label 'cpu_small'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(barcodes), path(adts), path(matrix), path(translate)
    val images_pulled

    output:
    tuple val(sampleId), path("*unfiltered.mtx.gz"), path("*unfiltered.barcodes.tsv"), path("*unfiltered.genes.tsv"), emit: merged_counts
    tuple val(sampleId), path("merged.mtx"), path("merged.barcodes.txt"), path("merged.genes.txt"), emit: old_merged_counts

    script:
    """
    publicCytometryMergeBarcodes.py -b ${barcodes} -a ${adts} -m ${matrix} -t ${translate} -s ${sampleId}
    """
}
