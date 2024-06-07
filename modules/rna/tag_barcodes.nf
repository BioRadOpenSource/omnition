/*
Tagging BAM with barcode features
*/

params.options = [:]

process TAG_BARCODES {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_large'
    label 'memory_small'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*_Aligned.sortedByCoord.tagged.bam'),
        path('*_Aligned.sortedByCoord.tagged.bam.bai'), emit: bam
    tuple val(sampleId), path('*_all_beads.tsv'), emit: beads

    script:
    """
    mkdir -p tmp/
    rnaBamTagger.py -b ${bam} -c ${task.cpus} -cp 0 -up 1 -ct XC -ut XM -t tmp/ -o ./
    ARGS=\$(cat tmp/samtools_cat.txt)
    samtools cat \$ARGS | samtools view -b -1 -o ${sampleId}_Aligned.sortedByCoord.tagged.bam -
    mv all_beads.tsv ${sampleId}_all_beads.tsv
    rm -r tmp/
    samtools index -@ ${task.cpus} ${sampleId}_Aligned.sortedByCoord.tagged.bam
    """
}
