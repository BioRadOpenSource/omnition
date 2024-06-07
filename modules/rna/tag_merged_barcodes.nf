/*
Tagging BAM with barcode features
*/

params.options = [:]

process TAG_MERGED_BARCODES {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(index), path(lookup)
    val images_pulled

    output:
    tuple val(sampleId), path('*.final.bam'),
        path('*.final.bam.bai'), emit: bam

    script:
    """
    mkdir -p tmp/
    rnaBamTagger.py -b ${bam} -c ${task.cpus} -ct XC -t tmp/ -o ./ -l ${lookup}
    ARGS=\$(cat tmp/samtools_cat.txt)
    samtools cat \$ARGS | samtools view -b -1 -o ${sampleId}.final.bam -
    rm -r tmp/

    samtools index -@ ${task.cpus} ${sampleId}.final.bam
    """
}
