/*
Merge chromosome fragments files back together
*/

process FINAL_FRAG_MERGE {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/deconvolution",
        pattern: '*.{fragments.tsv.gz,fragments.tsv.gz.tbi}', mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bedpe)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.fragments.tsv.gz"),
    path("${sampleId}.fragments.tsv.gz.tbi"), emit: final_frag_merge

    script:
    """
    cat ${bedpe} | cut -f 1,2,3,4,6 > ${sampleId}.fragments.tsv
    bgzip ${sampleId}.fragments.tsv
    tabix -p bed ${sampleId}.fragments.tsv.gz
    """
}
