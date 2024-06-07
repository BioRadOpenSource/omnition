/*
Aggregating all sample metrics together into a tidy format
*/

params.options = [:]

process AGGREGATE_METRICS {
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}/", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'
    if (!(workflow.profile =~ /(awsbatch|tower)/)) {
        beforeScript 'ulimit -Ss unlimited'
    }

    input:
    path input
    val assay
    val images_pulled

    output:
    path 'metric_summary.csv', emit: summary

    script:
    """
    atacAggregateMetrics.R ./ "${assay}" "${params.options.mitoContig}" "${params.catac != null}"
    """
}
