/*
Update single bead cell barcodes to proper format
*/

params.options = [:]

process CORRECT_CELL_BARCODES {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(barcode_translate), path(filtered_barcodes)
    val images_pulled

    output:
    tuple val(sampleId), path("*_barcodeTranslate.tsv"), emit: barcode_translate

    script:
    """
    rnaCorrectCellBarcodes.py -b ${filtered_barcodes} \
        -bt ${barcode_translate} \
        -s ${sampleId}
    """
}
