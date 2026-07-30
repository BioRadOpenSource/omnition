/*
Aggregating read counts per gene per cell into a count matrix per sample
*/

process MAKE_COUNT_MATRIX {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/countMatrix/",
        mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(count), path(symbols)
    val images_pulled

    output:
    tuple val(sampleId), path('*.unfiltered.mtx.gz'), path('*.unfiltered.barcodes.tsv'),
        path('*.unfiltered.genes.tsv'), emit: matrix

    script:
    """
    publicCoreBuildCountMatrix.R -o ./ -p ${sampleId}.unfiltered -w NULL -g ${symbols} ${count}
    pigz -p ${task.cpus} ${sampleId}.unfiltered.mtx
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=4

    # Check that we have expected number of sample input files
    check_input_files "MAKE_COUNT_MATRIX" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.unfiltered.mtx.gz" "${sampleId}.unfiltered.barcodes.tsv"
    "${sampleId}.unfiltered.genes.tsv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("count" "task.cpus" "workflow.manifest.version" "params.options.resultsDir"
    "sampleId" "symbols")

    parameters=("${count}" "${task.cpus}" "${workflow.manifest.version}"
    "${params.options.resultsDir}" "${sampleId}" "${symbols}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MAKE_COUNT_MATRIX" "\$param_name" "\$parameter"
    done
    """
}
