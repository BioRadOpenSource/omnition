/*
Count the frequency of different features for mixed species
*/

process MIXED_EXPRESSION {
    tag "${sampleId}-${species}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), val(species), path(allowlist), path(symbols), path(count)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.${species}_matrix_features.csv"), emit:mixed_features

    script:
    """
    publicCoreBuildCountMatrix.R -o ./ -p ${sampleId}.${species}.filtered -w ${allowlist} -g ${symbols} ${count} -m mixed
    pigz -p ${task.cpus} ${sampleId}.${species}.filtered.mtx

    publicRnaCountMatrixFeatures.R ${sampleId}.${species} ${sampleId}.${species}.filtered.mtx.gz
    """

    stub:
    """
    # Create expected output files
    outputs=("${sampleId}.${species}_matrix_features.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("task.cpus" "symbols" "sampleId" "workflow.manifest.version" "species"
    "allowlist" "count")

    parameters=("${task.cpus}" "${symbols}" "${sampleId}" "${workflow.manifest.version}"
    "${species}" "${allowlist}" "${count}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MIXED_EXPRESSION" "\$param_name" "\$parameter"
    done
    """
}
