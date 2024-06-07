/*
Make the count matrix
*/

params.options = [:]

process MAKE_COUNT_MATRIX {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/countMatrix",
        pattern: '*{.count_matrix.mtx.gz,.row_names.txt.gz,.column_names.txt.gz}', mode: 'copy',
        overwrite: true
    publishDir "${params.options.resultsDir}/${sampleId}/alignments",
        pattern: '*.cell_data.csv', mode: 'copy',
        overwrite: true, enabled: params.catac != null
    publishDir "${params.options.resultsDir}/${sampleId}/cellFilter",
        pattern: '*.cell_data.csv', mode: 'copy', // groovylint-disable-line
        overwrite: true, enabled: params.catac == null
    label 'cpu_xsmall'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(data), path(bam), path(bai), path(bed)
    val images_pulled

    output:
    tuple path("${sampleId}.column_names.txt.gz"), path("${sampleId}.row_names.txt.gz"),
        path("${sampleId}.count_matrix.mtx.gz"), path("${sampleId}.read_counts.txt"), emit: matrix
    tuple val(sampleId), path("${sampleId}.count_matrix.mtx.gz"), emit:mtx
    tuple val(sampleId), path("${sampleId}.cell_data.csv"), emit:stats

    script:
    mixedcheck = params.options.mixed ? "--mixed" : ""
    """
    atacCountMatrix.py -b ${bam} -w ${data} -p ${bed} -c ${task.cpus} -s ${sampleId} -o ./ ${mixedcheck}
    pigz -p ${task.cpus} ${sampleId}.count_matrix.mtx
    pigz -p ${task.cpus} ${sampleId}.row_names.txt
    pigz -p ${task.cpus} ${sampleId}.column_names.txt
    """
}
