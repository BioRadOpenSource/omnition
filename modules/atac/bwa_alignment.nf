/*
Align ATAC reads to reference genome
*/

params.options = [:]

process BWA_ALIGNMENT {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_xlarge'
    label 'memory_xlarge'

    input:
    tuple val(sampleId), path(fastq), path(reference)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.alignments.tagged.bam"),
        path("${sampleId}.alignments.tagged.bam.bai"), emit: bam
    tuple val(sampleId), path('*_bwa_read_counts.csv'), emit: count

    script:
    """
    bwa mem -M -t ${task.cpus} -L 5 bwa-index/\$(basename \$(ls bwa-index | head -n1) .amb) ${fastq[0]} ${fastq[1]} |
        atacBamTagger.py -cp 0 -ct XB |
        samtools sort -@ ${task.cpus} -O bam -o ${sampleId}.alignments.tagged.bam
    samtools index -@ ${task.cpus} ${sampleId}.alignments.tagged.bam

    # Calculate the number of output reads by dividing the number of mapped proper pairs by 2
    proper_pairs=\$(samtools view -c -f 3 -F 4 -F 256 ${sampleId}.alignments.tagged.bam)
    out_count=\$((\${proper_pairs} / 2))

    # Write number of input and output reads to file
    READCOUNTFILE="${sampleId}_bwa_read_counts.csv"
    printf "sample,process,metric,count\n" > \$READCOUNTFILE
    printf "${sampleId},bwa,input,\$(zcat ${fastq[0]} | awk 'END {print NR/4}')\n" >> \$READCOUNTFILE
    printf "${sampleId},bwa,output,\${out_count}\n" >> \$READCOUNTFILE
    """
}
