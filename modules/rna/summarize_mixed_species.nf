/*
Calculating crosstalk statistics, etc. for mixed species experiments
*/

params.options = [:]

process SUMMARIZE_MIXED_SPECIES {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/countMatrix/", pattern:'*.rds', mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(species1), val(species2)
    tuple val(sampleId), path(counts), path(allowlist)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.species_mix_counts.rds"), emit: count
    tuple val(species1), val(sampleId), path("${sampleId}.${species1}.allowlist.csv"), emit: s1_allowlist
    tuple val(species2), val(sampleId), path("${sampleId}.${species2}.allowlist.csv"), emit: s2_allowlist
    tuple val(sampleId), path("${sampleId}_crosstalk_density.csv"), emit: crosstalk_density

    script:
    """
    coreSpeciesMixStats.R -p ${sampleId} -a ${allowlist} -ct ${params.options.crosstalkThreshold} \
        ${counts} ${species1} ${species2}
    """
}
