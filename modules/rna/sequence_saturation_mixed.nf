/*
Calculate sequence saturation for each species in each sample
*/

params.options = [:]

process SEQUENCE_SATURATION_MIXED {
    tag "${sampleId}-${species}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), val(species), path(allowlist), path(symbols), path(fragments), path(summary)
    val assay
    val images_pulled

    output:
    tuple val(sampleId), path("*_sequence_saturation.csv"), emit: results

    script:
    assay = (assay == "3' RNA Droplet") ? "RNA" : assay
    """
    coreSequenceSaturation.R \
        --sample_id ${sampleId} \
        --feature_counts_file ${fragments} \
        --metrics_file ${summary} \
        --allowlist_file ${allowlist} \
        --symbols_file ${symbols} \
        --species_id ${species} \
        --assay "${assay}"
    """
}
