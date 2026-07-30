/*
Clean summit peaks
*/

process CLEAN_PEAKS {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/peaks",
        pattern: '*.fixedwidthpeaks.bed', mode: 'copy', overwrite: true
    errorStrategy 'ignore'
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(bed), path(blocklist), path(sizes)
    val images_pulled

    output:
    tuple val(sampleId), path('*.fixedwidthpeaks.bed'), emit: peaks // groovylint-disable-line
    tuple val(sampleId), path('*.peak_calling_summary.csv'), emit: summary // groovylint-disable-line

    script:
    """
    publicAtacSummitsToCleanPeaks.R --peak_width 250 \
        --fdr_threshold 0.01 \
        --output_directory ./ \
        --name ${sampleId} \
        --blocklist ${blocklist} \
        ${bed} \
        ${sizes}
    """
}
