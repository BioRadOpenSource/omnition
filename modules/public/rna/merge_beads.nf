/*
Parsing and correcting cell barcodes
*/

process MERGE_BEADS {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(edges), path(beads)
    val images_pulled

    output:
    tuple val(sampleId), path("*.barcodeTranslate.tsv"), emit: barcode_translate
    tuple val(sampleId), path("*filtered_bead_list.csv"), emit: filtered_beads

    script:
    """
    # create barcode translation table
    publicRnaMergeBeads.py -e ${edges} \
        -s ${sampleId} \
        -u ${params.options.beadMergeUmiThreshold.get( sampleId )} \
        -bf ${beads} \
        -v
    mv filtered_bead_list.csv ${sampleId}.filtered_bead_list.csv
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=4

    # Check that we have expected number of sample input files
    check_input_files "MERGE_BEADS" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.barcodeTranslate.tsv" "${sampleId}.filtered_bead_list.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("beads" "oneToOne" "edges" "params.options.resultsDir" "sampleId"
    "workflow.manifest.version" "params.options.beadMergeUmiThreshold.get( sampleId )")

    parameters=("${beads}" "${oneToOne}" "${edges}" "${params.options.resultsDir}"
    "${sampleId}" "${workflow.manifest.version}"
    "${params.options.beadMergeUmiThreshold.get( sampleId )}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MERGE_BEADS" "\$param_name" "\$parameter"
    done
    """
}
