/*
Parsing and correcting cell barcodes
*/

process DEBARCODER {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
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
    r2_json = "/opt/biorad/assets/dead/${params.options.bead.format.get(sampleId)}".replaceAll(/r1$/, "r2.json")
    if (params.rna) {
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
    } else if (params.atac) {
        """
        export RUST_LOG=info
        dead ${fastq[0]} ${fastq[1]} \
            -c /opt/biorad/assets/dead/${params.options.bead.format.get( sampleId )}.json \
            -o ./ \
            -i ${sampleId}
        """
    } else {
        """
        # Parsing and correcting bead oligo
        # NOTE: '-c' path refers to path inside container

        dead ${fastq} -c /opt/biorad/assets/dead/shasta_normalbead.json \
            -a Edge \
            -s /opt/biorad/assets/dead/shasta_normalbead.json \
            -o ./ \
            -i ${sampleId} \
            --output_format fastq
        mv edges.tsv ${sampleId}_edges.tsv
        """
    }

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=7

    # Check that we have expected number of sample input files
    check_input_files "DEBARCODER" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_R1_barcode_stats.tsv" "${sampleId}_R2_barcode_stats.tsv"
    "${sampleId}_R1.trimmed_debarcoded.fastq.gz"
    "${sampleId}_R2.trimmed_debarcoded.fastq.gz" "${sampleId}_edges.tsv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.resultsDir" "workflow.manifest.version"
    "sampleId" "fastq" "params.options.bead.format.get( sampleId )")

    parameters=("${params.options.resultsDir}"
    "${workflow.manifest.version}" "${sampleId}" "${fastq}"
    "${params.options.bead.format.get( sampleId )}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "DEBARCODER" "\$param_name" "\$parameter"
    done
    """
}

