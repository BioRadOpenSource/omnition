/*
Merge chromosome bam files back together
*/

params.options = [:]

process FINAL_BAM_MERGE {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/alignments", pattern:'*.{bam,bai}',
        mode: 'copy', enabled: params.catac == null, overwrite: true
    label 'cpu_xlarge'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.final.bam"), path("${sampleId}.final.bam.bai"), emit: final_bam
    tuple val(sampleId), path('*_deconvolution_output_read_counts.csv'), emit: count

    script:
    """
    (samtools merge -@ ${task.cpus} ${sampleId}.final.bam ${bam}) &> ${sampleId}_finalmerge.log
    samtools index -@ ${task.cpus} ${sampleId}.final.bam

    # Count input Proper Pairs
    input_proper_pairs=0
    for FILE in ${bam}; do
        input_proper_pairs=\$((\$(samtools flagstat -@ ${task.cpus} \$FILE | grep "properly paired" \
        | awk '{print \$1}')+\${input_proper_pairs}))
    done
    final_in_count=\$((\${input_proper_pairs} / 2))

    # Count output proper pairs
    output_proper_pairs=\$(samtools flagstat -@ ${task.cpus} ${sampleId}.final.bam \
    | grep "properly paired" | awk '{print \$1}')
    out_count=\$((\${output_proper_pairs} / 2))

    # Write number of output reads to file
    READCOUNTFILE="${sampleId}_deconvolution_output_read_counts.csv"
    printf "sample,process,metric,count\n" > \$READCOUNTFILE
    printf "${sampleId},deconvolution,input,\${final_in_count}\n" >> \$READCOUNTFILE
    printf "${sampleId},deconvolution,output,\${out_count}\n" >> \$READCOUNTFILE
    """
}
