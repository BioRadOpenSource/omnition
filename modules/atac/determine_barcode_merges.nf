/*
Determining which barcodes to merge
*/

params.options = [:]

process DETERMINE_BARCODE_MERGES {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"

    if (workflow.profile == 'aws') {
        label 'xlarge'
  } else {
        label 'cpu_medium'
        label 'memory_small'
    }

    publishDir "${params.resultsDir}/${sampleId}/deconvolution", pattern: '*.barcodeTranslate.tsv',
      mode: 'copy', enabled: !params.options.barcodedTn5, overwrite: true

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
  MINIMUM_JACCARD_INDEX=0.0
  ONE_TO_ONE=False
  BARCODE_PRIOR_FILE="none"
  if [ "${params.options.barcodedTn5}" = "true" ]; then
    BARCODED_TN5=True

    # merge together fastq-ti files to get fastq files
    cat ${quant} > merged.barcodeQuantSimple.csv
    cat ${beads} > merged_barcode_allowlist.csv

    atacCallDoublets.R \
      ./ \
      merged.barcodeQuantSimple.csv \
      merged_barcode_allowlist.csv \
      new.deconvolutionParams.tmp \
      ${sampleId}.implicatedBarcodes.csv.gz \
      \$MINIMUM_JACCARD_INDEX \
      ${sampleId} \
      \$ONE_TO_ONE \
      \$BARCODED_TN5 \
      ${tilen} \
      \$BARCODE_PRIOR_FILE

    # make the bapParams files
    for FILE in ${param}; do
        # cutting out the ".orig" part of the file name of the input file
        params_file=\${FILE:0:-9}\${FILE: -4}
        cat \$FILE new.deconvolutionParams.tmp > \$params_file
    done

  else
    BARCODED_TN5=False

    atacCallDoublets.R \
      ./ \
      ${quant} \
      ${beads} \
      new.deconvolutionParams.tmp \
      ${sampleId}.implicatedBarcodes.csv.gz \
      \$MINIMUM_JACCARD_INDEX \
      ${sampleId} \
      \$ONE_TO_ONE \
      \$BARCODED_TN5 \
      ${tilen} \
      \$BARCODE_PRIOR_FILE
    cat ${param} new.deconvolutionParams.tmp > "${sampleId}.deconvolutionParams.csv"
  fi

  """
}
