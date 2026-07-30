/*
Calculate sequence saturation for each species in each sample
*/

process SEQUENCE_SATURATION_MIXED {
    tag "${sampleId}-${species}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), val(species), path(allowlist), path(symbols), path(fragments), path(summary)
    val assay
    val images_pulled

    output:
    tuple val(sampleId), path("*_sequence_saturation.csv"), emit: results

    script:
    """
    publicCoreSequenceSaturation.R \
        --sample_id ${sampleId} \
        --feature_counts_file ${fragments} \
        --metrics_file ${summary} \
        --allowlist_file ${allowlist} \
        --symbols_file ${symbols} \
        --species_id ${species} \
        --assay "${assay}"
    """

    stub:
    """
    # Create expected output files
    outputs=("${sampleId}.${species}_sequence_saturation.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("sampleId" "symbols" "workflow.manifest.version" "species" "allowlist"
    "summary" "assay" "cDNA_counts" "fragments")

    parameters=("${sampleId}" "${symbols}" "${workflow.manifest.version}" "${species}"
    "${allowlist}" "${summary}" "${assay}" "${cDNA_counts}" "${fragments}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "SEQUENCE_SATURATION_MIXED" "\$param_name" "\$parameter"
    done
    """
}
