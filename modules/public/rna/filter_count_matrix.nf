/*
Filter count matrix to include only barcodes on the allowlist
*/

process FILTER_COUNT_MATRIX {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/countMatrix/",
        mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(count), path(allowlist), path(symbols)
    val images_pulled

    output:
    tuple val(sampleId), path('*.filtered.mtx.gz'), path('*.filtered.barcodes.tsv'),
        path('*.filtered.genes.tsv'), emit: matrix

    script:
    """
    publicCoreBuildCountMatrix.R -o ./ -p ${sampleId}.filtered -w ${allowlist} -g ${symbols} ${count}
    pigz -p ${task.cpus} ${sampleId}.filtered.mtx
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=5

    # Check that we have expected number of sample input files
    check_input_files "FILTER_COUNT_MATRIX" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.filtered.mtx.gz" "${sampleId}.filtered.barcodes.tsv"
    "${sampleId}.filtered.genes.tsv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("task.cpus" "count" "params.options.resultsDir" "workflow.manifest.version"
    "allowlist" "symbols" "sampleId")

    parameters=("${task.cpus}" "${count}" "${params.options.resultsDir}"
    "${workflow.manifest.version}" "${allowlist}" "${symbols}" "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "FILTER_COUNT_MATRIX" "\$param_name" "\$parameter"
    done
    """
}
