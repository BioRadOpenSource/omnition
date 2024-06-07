/*
Tagging BAM with genomic features
*/

params.options = [:]

process TAG_FEATURES {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.featuresTagged.bam'), path('*.featuresTagged.bam.bai'), emit: bam

    script:
    include_introns = params.options.includeIntrons == true ? "--include-introns" : ""
    """
    mkdir -p tmp/
    rnaGeneTagger.py -b ${bam} -gt XT -ft XF -t tmp/ -o ./ ${include_introns}
    rm -r tmp
    mv tagged.bam ${sampleId}.featuresTagged.bam

    samtools index -@ ${task.cpus} ${sampleId}.featuresTagged.bam
    """
}
