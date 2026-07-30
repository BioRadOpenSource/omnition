/*
Tagging BAM with barcode features
*/

process TAG_MERGED_BARCODES {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/alignment/",
        pattern:'*{.final.bam,.final.bam.bai}', mode: 'copy', overwrite: true
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(index), path(lookup)
    val images_pulled

    output:
    tuple val(sampleId), path('*.final.bam'),
        path('*.final.bam.bai'), emit: bam

    script:
    """
    mkdir -p tmp/
    omnition rna tag-bam \
        --bam-file ${bam} --n-threads ${task.cpus} --cell-barcode-tag XC \
        --cell-barcode-position 0 \
        --umi-position 1 --umi-tag XM --tmp-dir tmp --lookup ${lookup}
    ARGS=\$(cat tmp/samtools_cat.txt)
    samtools cat \$ARGS | samtools view -b -1 -o ${sampleId}.final.bam -
    rm -r tmp/
    # removing file that gets created for tag_barcodes
    rm all_beads.tsv

    samtools index -@ ${task.cpus} ${sampleId}.final.bam
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=4

    # Check that we have expected number of sample input files
    check_input_files "TAG_MERGED_BARCODES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.final.bam" "${sampleId}.final.bam.bai"
    "${sampleId}.final.bam.bai")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("lookup" "workflow.manifest.version" "bam" "task.cpus" "sampleId")

    parameters=("${lookup}" "${workflow.manifest.version}" "${bam}" "${task.cpus}"
    "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "TAG_MERGED_BARCODES" "\$param_name" "\$parameter"
    done
    """
}
