/*
Aggregating all sample metrics together into a tidy format
*/

params.options = [:]

process AGGREGATE_METRICS {
    beforeScript 'ulimit -Ss unlimited'
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.reportsDir}/", mode: 'copy', overwrite: true
    if (workflow.profile == 'aws') {
        label 'medium'
    } else {
        label 'cpu_small'
        label 'memory_xsmall'
    }

    input:
    path input
    val assay
    val images_pulled

    output:
    path 'metric_summary.csv', emit: summary

    script:
    """
    atacAggregateMetrics.R ./ "${assay}" "${params.atac.mitoContig}" "${params.options.barcodedTn5}"
    """
}
