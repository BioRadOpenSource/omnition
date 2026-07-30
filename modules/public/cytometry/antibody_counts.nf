/*
Counting the antibodies per cell in the fully processed cells
*/

process ANTIBODY_COUNTS {
    tag "${sampleId}"
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/averageCounts/",
        mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(fully_processed)
    val images_pulled

    output:
    path "*_antibody_counts_per_cell.csv", emit: antibody_counts

    script:
    """
    publicCytometryAntibodiesPerCell.py -c ${fully_processed} -s ${sampleId}
    """
}
