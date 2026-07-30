/*
Generate pipeline summary
*/

process BUILD_REPORT_CONTENTS {
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}", mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_xxsmall'

    input:
    path input
    val images_pulled

    output:
    path("*_summary_table.csv"), emit: pipeline_summary
    path("fastqTIreadcountsfinal.csv"), emit: index_counts optional true

    script:
    """
    # If a file starting with metric_summary_ exists, copy its contents to metric_summary.csv
    if ls metric_summary_* 1> /dev/null 2>&1; then
        cat metric_summary_* > metric_summary.csv
    fi

    mkdir tmp
    echo ${params} | sed 's/\\[/\\[\\n/g' | sed 's/,/\\n/g' | sed 's/:/,/g' |
        tr -d '[' | tr -d ']' | tail -n +3 > ./tmp/params.csv

    publicAtacBuildReportContents.R ./ ./ "${params.catac != null && params.atac == null}"
    """
}
