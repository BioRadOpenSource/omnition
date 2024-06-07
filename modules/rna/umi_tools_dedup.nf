/*
Removing duplicate reads
*/

params.options = [:]

process UMI_TOOLS_DEDUP {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/alignment/",pattern:'*.{bam,bai}',
        mode: 'copy', overwrite: true
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
    PYTHONHASHSEED=0 umi_tools dedup \
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
}
