/*
Annotate barcodes, need to add script still.
*/

params.options = [:]

process TAG_BARCODES {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'xlarge'
    } else {
        label 'cpu_xlarge'
        label 'memory_xlarge'
    }

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.alignments.tagged.bam"),
        path("${sampleId}.alignments.tagged.bam.bai"), emit: bam

    script:
    """
    mkdir -p tmp/
    atacBamTagger.py -b ${bam} -c ${task.cpus} -cp 0 -ct XB -t tmp/ -o ./
    rm -r tmp/
    mv tagged.bam ${sampleId}.alignments.tagged.tmp.bam
    samtools sort -@ ${task.cpus} ${sampleId}.alignments.tagged.tmp.bam -o ${sampleId}.alignments.tagged.bam
    samtools index -@ ${task.cpus} ${sampleId}.alignments.tagged.bam
    """
}
