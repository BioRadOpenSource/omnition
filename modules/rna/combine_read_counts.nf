/*
Combine read counts files
*/

params.options = [:]

process COMBINE_BEAD_COUNTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(counts_duplicated), path(counts_unmapped, stageAs: 'unmapped.sequence_counts.csv')
    val images_pulled

    output:
    tuple val(sampleId), path('*.counts_per_barcode.csv'), emit: sequence_counts

    script:
    """
    # Generate a per-bead per-contig count table
    rnaCombineReadCounts.R ${counts_duplicated} unmapped.sequence_counts.csv -p ${sampleId}
    """
}
