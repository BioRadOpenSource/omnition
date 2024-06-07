/*
Merge read counts from reannotate bam back together
*/

params.options = [:]

process MERGE_REANN_READ_COUNTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(counts_files)
    val images_pulled

    output:
    tuple val(sampleId), path('*_merged_reannotate_read_counts.csv'), emit: count

    script:
    '''
    atacMergeReannCounts.R ./
    '''
}
