/*
Make the count matrix
*/

params.options = [:]

process MAKE_COUNT_MATRIX {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'large'
    } else {
        label 'cpu_xsmall'
        label 'memory_xlarge'
    }

    publishDir "${params.resultsDir}/${sampleId}/countMatrix",
        pattern: '*{.count_matrix.mtx.gz,.row_names.txt.gz,.column_names.txt.gz}', mode: 'copy',
        overwrite: true

    input:
    tuple val(sampleId), path(data), path(bam), path(bai), path(bed)
    val images_pulled

    output:
    tuple path("${sampleId}.column_names.txt.gz"), path("${sampleId}.row_names.txt.gz"),
        path("${sampleId}.count_matrix.mtx.gz"), path("${sampleId}.read_counts.txt"), emit: matrix
    tuple val(sampleId), path("${sampleId}.count_matrix.mtx.gz"), emit:mtx

    script:
    """
    atacCountMatrix.py -b ${bam} -w ${data} -p ${bed} -c ${task.cpus} -s ${sampleId} -o ./
    pigz -p ${task.cpus} ${sampleId}.count_matrix.mtx
    pigz -p ${task.cpus} ${sampleId}.row_names.txt
    pigz -p ${task.cpus} ${sampleId}.column_names.txt
    """
}
