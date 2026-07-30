/*
Calculate mixed species metrics
*/

process SUMMARIZE_MIXED_SPECIES {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/alignments",
        pattern: '*tsv', mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(data)
    tuple env(species1), env(species2)
    val images_pulled

    output:
    tuple val(sampleId), path('*.crosstalk.csv'), path('*.species_mix_counts.csv'), emit: stats

    script:
    """
    publicCoreSpeciesMixStats.R -p ${sampleId} -y "atac" -ct ${params.options.crosstalkThreshold} \
        ${data} \$species1 \$species2
    """
}
