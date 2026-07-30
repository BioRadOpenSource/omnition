/*
Removing duplicate reads
*/

process UMI_TOOLS_DEDUP {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.deduped.bam'), path('*.deduped.bam.bai'), emit: bam
    tuple val(sampleId), path('*_umi_tools_dedup_read_counts.csv'), emit: count

    script:
    """
    samtools sort -n -@ ${task.cpus} ${bam} | \
    samtools sort -@ ${task.cpus} - | \
    samtools view --write-index -@ ${task.cpus} -F 256 -o ${sampleId}.nosecondary.bam -

    mkdir -p tmp/

    # Deduplicating bam file
    umi_tools dedup \
        --temp-dir "tmp/" \
        --per-cell \
        --per-gene \
        --extract-umi-method tag \
        --cell-tag XC \
        --umi-tag XM \
        --gene-tag XT \
        --random-seed 1 \
        --skip-tags-regex "^ " \
        -I ${sampleId}.nosecondary.bam \
        -L ${sampleId}.dedup.log \
        -S ${sampleId}.deduped.bam \
        --compress 2

    # Indexing deduped bam file
    samtools index -@ ${task.cpus} ${sampleId}.deduped.bam

    # Write number of input and output reads to file
    READCOUNTFILE="${sampleId}_umi_tools_dedup_read_counts.csv"
    printf "sample,process,metric,value\n" > \$READCOUNTFILE
    printf "${sampleId},umi_tools_dedup,input,\$(samtools view -F 256 -c ${bam})\n" >> \$READCOUNTFILE
    printf "${sampleId},umi_tools_dedup,output,\$(samtools view -F 256 -c \
        ${sampleId}.deduped.bam)\n" >> \$READCOUNTFILE
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=8

    # Check that we have expected number of sample input files
    check_input_files "UMI_TOOLS_DEDUP" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.deduped.bam" "${sampleId}.deduped.bam.bai"
    "${sampleId}.deduped.bam.bai" "${sampleId}_umi_tools_dedup_read_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.resultsDir" "bam" "task.cpus" "workflow.manifest.version"
    "sampleId")

    parameters=("${params.options.resultsDir}" "${bam}" "${task.cpus}"
    "${workflow.manifest.version}" "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "UMI_TOOLS_DEDUP" "\$param_name" "\$parameter"
    done
    """
}
