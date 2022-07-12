/*
Generate output reports
*/

params.options = [:]

process BUILD_REPORT_CONTENTS {
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'small'
    } else {
        label 'cpu_medium'
        label 'memory_xxsmall'
    }

    publishDir "${params.reportsDir}", mode: 'copy'

    input:
    path input
    val images_pulled

    output:
    path("*_summary_table.csv"), emit: pipeline_summary
    path("metric_summary_updated.csv"), emit: metric_summary
    path("fastqTIreadcountsfinal.csv"), emit: index_counts optional true

    script:
    """
    mkdir tmp
    echo ${params} | sed 's/,/\\n/g' | sed 's/:/,/g' | tr -d '[' | tr -d ']' | tail -n +3 > ./tmp/params.csv

    atacBuildReportContents.R ./ ./ ${params.options.barcodedTn5}
    """
}
