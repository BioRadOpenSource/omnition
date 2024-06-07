/*
Parsing and correcting cell barcodes
*/

params.options = [:]

process DEBARCODER {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(fastq)
    val images_pulled

    output:
    tuple val(sampleId), path('*_R1_barcode_stats.tsv'), emit: count
    tuple val(sampleId), path('*_R2_barcode_stats.tsv'), emit: r2_stats
    tuple val(sampleId), path('*_debarcoded.fastq.gz'), emit: fastq
    tuple val(sampleId), path('*_edges.tsv'), optional: true, emit: edges

    script:
    """
    # Parsing and correcting bead oligo
    # NOTE: '-c' path refers to path inside container

    dead ${fastq} -c /opt/biorad/assets/dead/${params.options.bead.format.get( sampleId )}.json \
        -a Edge \
        -s /opt/biorad/assets/dead/${params.options.bead.format.get( sampleId )}.json \
        -o ./ \
        -i ${sampleId}
    mv edges.tsv ${sampleId}_edges.tsv
    """
}
