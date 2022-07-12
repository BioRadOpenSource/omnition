/*
Calculate mixed species metrics
*/

params.options = [:]

process SUMMARIZE_MIXED_SPECIES {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'small'
    } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }

    input:
    tuple val(sampleId), path(data)
    tuple env(species1), env(species2)
    val images_pulled

    output:
    tuple val(sampleId), path('*.crosstalk.csv'), path('*.species_mix_counts.csv'), emit: stats

    script:
    """
    coreSpeciesMixStats.R -p ${sampleId} -y "atac" -ct ${params.options.crosstalkthreshold} \
        ${data} \$species1 \$species2
    """
}
