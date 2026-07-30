/*
Aligning reads to the reference genome(s)
*/

process STAR_ALIGN {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xlarge'
    label 'memory_large'

    input:
    tuple val(sampleId), path(fastq), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*_Aligned.out.bam'), emit: raw_bam
    tuple val(sampleId), path('*Log.final.out'), emit: log

    script:
    """
    STAR \
        --runThreadN ${task.cpus} \
        --limitIObufferSize 60000000 \
        --limitOutSJcollapsed 2000000 \
        --genomeDir ${index} \
        --readFilesIn ${fastq[1]} \
        --readFilesCommand zcat \
        --outFileNamePrefix ${sampleId}_ \
        --outSAMunmapped Within \
        --outSAMtype BAM Unsorted \
        --outBAMcompression 2
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=7

    # Check that we have expected number of sample input files
    check_input_files "STAR_ALIGN" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_Aligned.out.bam" "${sampleId}_Log.final.out")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "task.cpus" "fastq[1]" "index" "sampleId")

    parameters=("${workflow.manifest.version}" "${task.cpus}" "${fastq[1]}" "${index}"
    "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "STAR_ALIGN" "\$param_name" "\$parameter"
    done
    """
}
