/*
Split Bam into separate species files
*/

params.options = [:]

process SPLIT_MIXED_BAM {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
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
    coreSplitMixedBam.py --id ${sampleId} --bam ${bam} -s1 ${species1} -s2 ${species2}

    # Write number of input and output reads to file
    READCOUNTFILE="${sampleId}_split_mixed_bam_star_align_read_counts.csv"
    printf "sample,process,read,metric,value\n" > \$READCOUNTFILE
    printf ""${sampleId}_${species1}",star_align,NA,filtered_reads,\$(grep 'input' \
    ${sampleId}_star_align_read_counts.csv | awk -F, '{print \$4}')\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species1}",star_align,NA,uniquely_mapped_reads,\$(samtools view -1 -c -q 255 \
        ${species1}.mixed.filtered.bam)\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species1}",star_align,NA,reads_mapped_to_multiple_loci,\$(samtools view -1 -F 256 \
        ${species1}.mixed.filtered.bam | grep -v -E -w 'NH:i:0|NH:i:1' \
        | wc -l)\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species2}",star_align,NA,filtered_reads,\$(grep 'input' \
    ${sampleId}_star_align_read_counts.csv | awk -F, '{print \$4}')\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species2}",star_align,NA,uniquely_mapped_reads,\$(samtools view -1 -c -q 255 \
        ${species2}.mixed.filtered.bam)\n" >> \$READCOUNTFILE
    printf ""${sampleId}_${species2}",star_align,NA,reads_mapped_to_multiple_loci,\$(samtools view -1 -F 256 \
        ${species2}.mixed.filtered.bam | grep -v -E -w 'NH:i:0|NH:i:1' \
        | wc -l)\n" >> \$READCOUNTFILE
    """
}
