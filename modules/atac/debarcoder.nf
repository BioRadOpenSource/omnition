/*
Parsing and correcting cell barcodes
*/

params.options = [:]

process DEBARCODER {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_small'

    input:
    tuple val(sampleId), path(fastq), path(config), path(fastq_valid)
    val images_pulled

    output:
    tuple val(sampleId), path('*_R1_barcode_stats.tsv'), emit: count
    tuple val(sampleId), path('*_debarcoded.fastq.gz'), emit: fastq

    script:
    if (params.catac && params.options.i7AsTi) {
        """
        export RUST_LOG=info
        dead -a DefaultParser \
            -c /opt/biorad/assets/dead/atac.json \
            -o ./ \
            --i7-as-ti \
            -x ${config} \
            -i ${sampleId} \
            ${fastq}
        """
    }
    else if (params.catac && !params.options.i7AsTi && params.options.tiRead == 'r2') {
        """
        export RUST_LOG=info
        dead -a PairedParser \
            -c /opt/biorad/assets/dead/atac.json \
            -s ${config} \
            -o ./ \
            -i ${sampleId} \
            ${fastq}
        """
    }
    else if (params.catac && !params.options.i7AsTi && params.options.tiRead == 'r1') {
        """
        export RUST_LOG=info
        dead -a DefaultParser \
            -c ${config} \
            -o ./ \
            -i ${sampleId} \
            ${fastq}
        """
    } else {
        """
        export RUST_LOG=info
        dead ${fastq[0]} ${fastq[1]} \
            -c /opt/biorad/assets/dead/atac.json \
            -o ./ \
            -i ${sampleId}
        """
    }
}
