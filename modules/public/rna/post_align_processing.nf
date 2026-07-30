/*
Aligning reads to the reference genome(s)
*/

process POST_ALIGN_PROCESSING {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_large'
    label 'memory_small'

    input:
    tuple val(sampleId), path(fastq), path(raw_bam), path(log)
    val images_pulled

    output:
    tuple val(sampleId), path('*_Aligned.sortedByCoord.mapped.out.bam'),
        path('*_Aligned.sortedByCoord.mapped.out.bam.bai'), emit: bam
    tuple val(sampleId), path('*_read_counts.csv'), emit: count
    tuple val(sampleId), path('*_unmapped.out.fastq'), emit: unmapped_fastq
    tuple val(sampleId), path('*_unmapped.out.bam'), path('*_unmapped.out.bam.bai'), emit: unmapped

    script:
    """
    samtools sort -@ ${task.cpus} ${sampleId}_Aligned.out.bam -o ${sampleId}_Aligned.sortedByCoord.out.bam

    samtools view -@ ${task.cpus} -b -f 4 ${sampleId}_Aligned.sortedByCoord.out.bam > ${sampleId}_unmapped.out.bam
    samtools index -@ ${task.cpus} ${sampleId}_unmapped.out.bam

    samtools bam2fq ${sampleId}_unmapped.out.bam > ${sampleId}_unmapped.out.fastq

    # Create the mapped BAM
    samtools view -@ ${task.cpus} -b -F 4 ${sampleId}_Aligned.sortedByCoord.out.bam \
    > ${sampleId}_Aligned.sortedByCoord.mapped.out.bam
    samtools index -@ ${task.cpus} ${sampleId}_Aligned.sortedByCoord.mapped.out.bam

    # Write number of input and output reads to file
    READCOUNTFILE="${sampleId}_star_align_read_counts.csv"
    INPUT=\$(zcat ${fastq[1]} | awk 'END {print NR/4}')
    MULTIMAP=\$(grep "Number of reads mapped to multiple loci" ${sampleId}_Log.final.out | awk '{print \$NF}')
    UNIQUEMAP=\$(grep "Uniquely mapped reads number" ${sampleId}_Log.final.out | awk '{print \$NF}')
    READSUM=\$((\$MULTIMAP+\$UNIQUEMAP))
    READMAP=\$(jq -n \$((\$MULTIMAP+\$UNIQUEMAP))/\$INPUT*100)
    READMAPCONF=\$(jq -n \$UNIQUEMAP/\$INPUT*100)
    printf "sample,process,metric,value\n" > \$READCOUNTFILE
    printf "${sampleId},star_align,input,\$INPUT\n" >> \$READCOUNTFILE
    printf "${sampleId},star_align,output,\$(samtools view -F 256 -c \
        ${sampleId}_Aligned.sortedByCoord.mapped.out.bam)\n" >> \$READCOUNTFILE
    printf "${sampleId},star_align,reads_mapped_genome,\$READMAP\n" >> \$READCOUNTFILE
    printf "${sampleId},star_align,reads_mapped_genome_confidently,\$READMAPCONF\n" >> \$READCOUNTFILE
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=13

    # Check that we have expected number of sample input files
    check_input_files "POST_ALIGN_PROCESSING" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_Aligned.sortedByCoord.mapped.out.bam"
    "${sampleId}_Aligned.sortedByCoord.mapped.out.bam.bai"
    "${sampleId}_Aligned.sortedByCoord.mapped.out.bam.bai"
    "${sampleId}_star_align_read_counts.csv" "${sampleId}_unmapped.out.fastq"
    "${sampleId}_unmapped.out.bam" "${sampleId}_unmapped.out.bam.bai"

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("task.cpus" "workflow.manifest.version" "sampleId"
    "params.options.resultsDir" "fastq[1]")

    parameters=("${task.cpus}" "${workflow.manifest.version}" "${sampleId}"
    "${params.options.resultsDir}" "${fastq[1]}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "POST_ALIGN_PROCESSING" "\$param_name" "\$parameter"
    done
    """
}
