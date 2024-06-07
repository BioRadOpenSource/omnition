/*
Determining which barcodes to merge
*/

params.options = [:]

process DETERMINE_BARCODE_MERGES {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/deconvolution", pattern: '*.barcodeTranslate.tsv',
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

    script:
    """
    if [ "${params.catac != null}" = "true" ]; then
        catac_assay=True

        # merge together fastq-ti files to get fastq files
        cat ${quant} > merged.barcodeQuantSimple.csv
        cat ${beads} > merged_barcode_allowlist.csv

        atacCallDoublets.py \
            -d ./ \
            -c merged.barcodeQuantSimple.csv \
            -q merged_barcode_allowlist.csv \
            -p new.deconvolutionParams.tmp \
            -n ${sampleId} \
            -b \$catac_assay \
            -l ${tilen}

        # make the bapParams files
        for FILE in ${param}; do
            # cutting out the ".orig" part of the file name of the input file
            params_file=\${FILE:0:-9}\${FILE: -4}
            cat \$FILE new.deconvolutionParams.tmp > \$params_file
        done
    else
        catac_assay=False

        atacCallDoublets.py \
            -d ./ \
            -c ${quant} \
            -q ${beads} \
            -p new.deconvolutionParams.tmp \
            -n ${sampleId} \
            -b \$catac_assay \
            -l ${tilen}
        cat ${param} new.deconvolutionParams.tmp > "${sampleId}.deconvolutionParams.csv"
    fi
    """
}
