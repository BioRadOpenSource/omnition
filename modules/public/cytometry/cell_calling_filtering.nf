/*
Cell Calling Filtering Process
*/

process CELL_CALLING_FILTERING {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/countMatrix/", mode: 'copy',
        overwrite: true, pattern: '*.filtered.*'
    label 'cpu_medium'
    label 'memory_large'

    input:
    tuple val(sampleId), path(allowlist), path(matrix), path(barcodes), path(genes)
    val images_pulled

    output:
    tuple val(sampleId), path("cell_calling.RDS"), path("umi_percell.RDS"), emit: cell_calling_filtering_results
    tuple val(sampleId), path("raw_cells.csv"), emit: raw_cells // groovylint-disable-line
    tuple val(sampleId), path("n_initial_cells.RDS"), emit: initial_cell_count
    tuple val(sampleId), path("deduped_reads.csv"), emit: deduped_reads
    tuple val(sampleId), path("filtered_cells.csv"), emit: filtered_cells
    tuple val(sampleId), path('*.filtered.mtx.gz'), path('*.filtered.barcodes.tsv'),
        path('*.filtered.genes.tsv'), emit: matrix
    tuple val(sampleId), path("*_cell_filtering.csv"), emit: filtering_metrics

    script:
    """
    publicCytometryCellFiltering.R \
        --inputdir \$PWD \
        --outputdir \$PWD \
        --sample ${sampleId} \
        --matrix_prefix "merged"
    pigz -p ${task.cpus} ${sampleId}.filtered.mtx
    """
}
