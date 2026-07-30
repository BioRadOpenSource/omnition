/*
Determining which barcodes to merge
*/

process DETERMINE_BARCODE_MERGES {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/deconvolution",
        pattern: '*.barcodeTranslate.tsv',
        mode: 'copy', enabled: params.catac == null, overwrite: true
    label 'cpu_large'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(quant), path(beads), path(param), path(overlap)
    val tilen
    val images_pulled

    output:
    path("*.implicatedBarcodes.csv.gz"), emit: implicated_barcodes
    path("*.barcodeTranslate.tsv"), emit: barcode_translate // groovylint-disable-line
    path("*.deconvolutionParams.csv"), emit: params
    path("*.deconvolution_metrics.csv"), emit: metrics

    script:
    catacBool = params.catac ? "--catac_assay" : ""
    """
    if [ "${params.catac != null}" = "true" ]; then
        catac_assay=True

        # merge together fastq-ti files to get fastq files
        cat ${quant} > merged.barcodeQuantSimple.csv
        cat ${beads} > merged_barcode_allowlist.csv
        publicAtacCallDoublets.py \
            -d ./ \
            -c merged.barcodeQuantSimple.csv \
            -q merged_barcode_allowlist.csv \
            -m 0 \
            -n ${sampleId} \
            ${catacBool} \
            -l ${tilen}
    else
        publicAtacCallDoublets.py \
            -d ./ \
            -c ${quant} \
            -q ${beads} \
            -m 0 \
            -n ${sampleId} \
            ${catacBool} \
            -l ${tilen}
    fi
    """
}
