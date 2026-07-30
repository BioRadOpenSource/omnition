/*
Combining library QC reports into a single report
*/

process COMPILE_QC_REPORT {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_small'

    input:
    path fastqc
    val images_pulled

    output:
    path 'multiqc_report.html', emit: qc_report

    script:
    """
    multiqc \
        --force \
        ./
    """
}
