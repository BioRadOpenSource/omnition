/*
Parsing and correcting cell barcodes
*/

params.options = [:]

process DEAD {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'medium'
    } else {
        label 'cpu_small'
        label 'memory_xxsmall'
    }

    input:
    tuple val(sampleId), path(fastq)
    path(config)
    val images_pulled

    output:
    tuple val(sampleId), path('*_barcode_stats.tsv'), emit: count
    tuple val(sampleId), path('*_debarcoded.fastq.gz'), emit: fastq

    script:
    if (params.options.barcodedTn5 && params.options.i7asti) {
        """
      export RUST_LOG=info
      dead -a DefaultParser -c /opt/biorad/assets/dead/atac.json -o ./ --i7-as-ti -x ${config} \
      -i ${sampleId} ${fastq}
      """
    }
    else if (params.options.barcodedTn5 && !params.options.i7asti && params.options.tiread == 'r2') {
        """
      export RUST_LOG=info
      dead -a PairedParser -c /opt/biorad/assets/dead/atac.json -s ${config} -o ./ -i ${sampleId} ${fastq}
      """
    }
    else if (params.options.barcodedTn5 && !params.options.i7asti && params.options.tiread == 'r1') {
        """
      export RUST_LOG=info
      dead -a DefaultParser -c ${config} -o ./ -i ${sampleId} ${fastq}
      """
    } else {
        """
      export RUST_LOG=info
      dead ${fastq[0]} ${fastq[1]} -c /opt/biorad/assets/dead/atac.json -o ./ -i ${sampleId}
      """
    }
}
