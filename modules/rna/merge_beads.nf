/*
Parsing and correcting cell barcodes
*/

params.options = [:]

process MERGE_BEADS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
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
    rnaMergeBeads.py -e ${edges} \
        -s ${sampleId} \
        -u ${params.options.beadMergeUmiThreshold.get( sampleId )} \
        -bf ${beads} \
        -v
    mv filtered_bead_list.csv ${sampleId}.filtered_bead_list.csv
    """
}
