/*
Merge chromosome fragments files back together
*/

params.options = [:]

process FINAL_FRAG_MERGE {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"

    if (workflow.profile == 'aws') {
        label 'small'
    } else {
        label 'cpu_medium'
        label 'memory_xsmall'
    }

    publishDir "${params.resultsDir}/${sampleId}/deconvolution",
      pattern: '*.{fragments.tsv.gz,fragments.tsv.gz.tbi}', mode: 'copy', overwrite: true

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
