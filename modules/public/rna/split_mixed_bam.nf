/*
Split Bam into separate species files
*/

process SPLIT_MIXED_BAM {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(bam_index), path(count)
    tuple val(species1), val(species2)
    val ch_images_pulled

    output:
    tuple val("${sampleId}_${species1}"), path("${species1}.mixed.filtered.bam"), emit: bam1
    tuple val("${sampleId}_${species2}"), path("${species2}.mixed.filtered.bam"), emit: bam2
    tuple val(sampleId), path('*_read_counts.csv'), emit: count

    script:
    """
    publicCoreSplitMixedBam.py --id ${sampleId} --bam ${bam} -s1 ${species1} -s2 ${species2}

    # Write number of input and output reads to file
    READCOUNTFILE="${sampleId}_split_mixed_bam_star_align_read_counts.csv"
    printf "sample,process,metric,value\n" > \$READCOUNTFILE
    printf ""${sampleId}_${species1}",star_align,filtered_reads,\$(grep 'input' \
    ${sampleId}_star_align_read_counts.csv | awk -F, '{print \$4}')\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species1}",star_align,uniquely_mapped_reads,\$(samtools view -1 -c -q 255 \
        ${species1}.mixed.filtered.bam)\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species1}",star_align,reads_mapped_to_multiple_loci,\$(samtools view -1 -F 256 \
        ${species1}.mixed.filtered.bam | grep -v -E -w 'NH:i:0|NH:i:1' \
        | wc -l)\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species2}",star_align,filtered_reads,\$(grep 'input' \
    ${sampleId}_star_align_read_counts.csv | awk -F, '{print \$4}')\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species2}",star_align,uniquely_mapped_reads,\$(samtools view -1 -c -q 255 \
        ${species2}.mixed.filtered.bam)\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species2}",star_align,reads_mapped_to_multiple_loci,\$(samtools view -1 -F 256 \
        ${species2}.mixed.filtered.bam | grep -v -E -w 'NH:i:0|NH:i:1' \
        | wc -l)\n" >> \$READCOUNTFILE
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=4

    # Check that we have expected number of sample input files
    check_input_files "SPLIT_MIXED_BAM" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${species1}.mixed.filtered.bam" "${species2}.mixed.filtered.bam"
    "${sampleId}_star_align_read_counts.csv"
    "${sampleId}_split_mixed_bam_star_align_read_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("species2" "bam" "sampleId" "species1" "workflow.manifest.version")

    parameters=("${species2}" "${bam}" "${sampleId}" "${species1}"
    "${workflow.manifest.version}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "SPLIT_MIXED_BAM" "\$param_name" "\$parameter"
    done
    """
}
