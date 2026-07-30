/*
Tagging BAM with barcode features
*/

process TAG_BARCODES {
    tag "${sampleId}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    label 'cpu_large'
    label 'memory_small'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*_Aligned.sortedByCoord.tagged.bam'),
        path('*_Aligned.sortedByCoord.tagged.bam.bai'), emit: bam

    script:
    umi_tag = params.options.assay == "imm_prof" ? "XU" : "XM"
    """
    mkdir -p tmp/
    omnition rna tag-bam --bam-file ${bam} \
        --n-threads ${task.cpus} \
        --cell-barcode-tag XC \
        --cell-barcode-position 0 \
        --umi-position 1 --umi-tag ${umi_tag} --tmp-dir tmp/
    ARGS=\$(cat tmp/samtools_cat.txt)
    samtools cat \$ARGS | samtools view -b -1 -o ${sampleId}_Aligned.sortedByCoord.tagged.bam -
    rm -r tmp/
    samtools index -@ ${task.cpus} ${sampleId}_Aligned.sortedByCoord.tagged.bam
    """

    stub:
    umi_tag = params.options.assay == "imm_prof" ? "XU" : "XM" //groovylint-disable-line
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=5

    # Check that we have expected number of sample input files
    check_input_files "TAG_BARCODES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_Aligned.sortedByCoord.tagged.bam"
    "${sampleId}_Aligned.sortedByCoord.tagged.bam.bai"
    "${sampleId}_Aligned.sortedByCoord.tagged.bam.bai" "${sampleId}_all_beads.tsv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("umi_tag" "task.cpus" "workflow.manifest.version" "bam" "sampleId")

    parameters=("${umi_tag}" "${task.cpus}" "${workflow.manifest.version}" "${bam}" "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "TAG_BARCODES" "\$param_name" "\$parameter"
    done
    """
}
