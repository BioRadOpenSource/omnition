/*
Update single bead cell barcodes to proper format
*/

process CORRECT_CELL_BARCODES {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(barcode_translate), path(filtered_barcodes)
    val images_pulled

    output:
    tuple val(sampleId), path("*_barcodeTranslate.tsv"), emit: barcode_translate

    script:
    """
    publicCoreCorrectCellBarcodes.py -b ${filtered_barcodes} \
        -bt ${barcode_translate} \
        -s ${sampleId}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=2

    # Check that we have expected number of sample input files
    check_input_files "CORRECT_CELL_BARCODES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_barcodeTranslate.tsv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "sampleId" "barcode_translate"
    "filtered_barcodes")

    parameters=("${workflow.manifest.version}" "${sampleId}" "${barcode_translate}"
    "${filtered_barcodes}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "CORRECT_CELL_BARCODES" "\$param_name" "\$parameter"
    done
    """
}
