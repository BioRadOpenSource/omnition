/*
Parsing and correcting cell barcodes
*/

process DEBARCODER_CITESEQ {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"

    label 'cpu_small'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(fastq), path(adt_config)
    val images_pulled

    output:
    tuple val(sampleId), path('*_R1_barcode_stats.tsv'), emit: count
    tuple val(sampleId), path('*_R2_barcode_stats.tsv'), emit: r2_stats
    tuple val(sampleId), path('*_R2_ADT_barcode_stats.tsv'), emit: adt_stats
    tuple val(sampleId), path('*_debarcoded.fastq.gz'), emit: fastq
    tuple val(sampleId), path('*_debarcoded_adt.fastq.gz'), emit: adt_fastq
    tuple val(sampleId), path('*_edges.tsv'), optional: true, emit: edges

    script:
    audit_flag = params.options.audit == true ? "-l" : ""
    """
    # Parsing and correcting bead oligo
    # NOTE: '-c' path refers to path inside container

    dead ${fastq} -c /opt/biorad/assets/dead/shasta_normalbead.json \
        -a EdgeADT \
        -s /opt/biorad/assets/dead/shasta_normalbead.json \
        -t ${adt_config}  \
        -o ./ \
        -i ${sampleId} \
        --output_format fastq
    mv edges.tsv ${sampleId}_edges.tsv
    """
}
