/*
Performs cell filtering and merging using RNA information for CiteSeq data.
*/

process CELL_FILTERING_MERGING {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/countMatrix/", mode: 'copy',
        overwrite: true, pattern: '*filtered.*'

    label 'cpu_medium'
    label 'memory_large'

    input:
    tuple val(sampleId), path(barcodes), path(adts), path(matrix), path(allowlist), path(translate)
    val images_pulled

    output:
    tuple val(sampleId), path("cell_calling.RDS"), path("umi_percell.RDS"), emit: cell_calling_filtering_results
    tuple val(sampleId), path("raw_cells.csv"), emit: raw_cells // groovylint-disable-line
    tuple val(sampleId), path("n_initial_cells.RDS"), emit: initial_cell_count
    tuple val(sampleId), path("deduped_reads.csv"), emit: deduped_reads
    tuple val(sampleId), path("filtered_cells.csv"), emit: filtered_cells
    tuple val(sampleId), path("*.filtered.mtx.gz"), path("*.filtered.barcodes.tsv"),
        path("*.filtered.genes.tsv"), emit: matrix
    tuple val(sampleId), path("*_cell_filtering.csv"), emit: filtering_metrics
    tuple val(sampleId), path("*.unfiltered.mtx.gz"), path("*.unfiltered.barcodes.tsv"),
        path("*.unfiltered.genes.tsv"), emit: unfiltered_matrix

    script:
    """
    publicCiteseqRnaCallingMerging.R \
        --sample_id ${sampleId} \
        --allowlist ${allowlist} \
        --translate ${translate}  \
        --matrix_prefix "merge"
    pigz -p ${task.cpus} ${sampleId}.filtered.mtx
    pigz -p ${task.cpus} ${sampleId}.unfiltered.mtx
    """
}
