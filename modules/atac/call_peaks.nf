/*
Call peaks
*/

params.options = [:]

process CALL_PEAKS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/peaks",
        pattern: '*_{summits.bed,peaks.narrowPeak,peaks.xls}', mode: 'copy',
        overwrite: true
    label 'cpu_xsmall'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    path('*')
    tuple val(sampleId), path('*_summits.bed'), emit:peaks

    script:
    """
    macs2 callpeak \
        --treatment ${bam} \
        --format BAMPE \
        --name ${sampleId} \
        --nomodel \
        --nolambda \
        --call-summits \
        --outdir ./
    """
}
