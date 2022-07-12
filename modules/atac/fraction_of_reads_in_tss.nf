/*
Calculate fraction of reads in TSS
*/

params.options = [:]

process FRACTION_OF_READS_IN_TSS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'small'
    } else {
        label 'cpu_xsmall'
        label 'memory_xsmall'
    }

    input:
    tuple val(sampleId), path(bam), path(index)
    path reference_tss
    val images_pulled

    output:
    tuple val(sampleId), path('*.frit.csv'), emit: frit

    script:
    """
    bedtools intersect -a ${bam} -b ${reference_tss} > readsInTss.bam
    readsInTSS=\$(samtools view -c readsInTss.bam)
    totalReads=\$(samtools view -c ${bam})
    tssReadPercentage=\$(printf "%.1f\n" \$((10**3 * \$readsInTSS/\$totalReads))e-1)
    printf "sample,process,metric,value\n" > ${sampleId}.frit.csv
    printf "${sampleId},frit,readsInTSS,\${readsInTSS}\n" >> ${sampleId}.frit.csv
    printf "${sampleId},frit,totalReads,\${totalReads}\n" >> ${sampleId}.frit.csv
    printf "${sampleId},frit,tssReadPercentage,\${tssReadPercentage}\n" >> ${sampleId}.frit.csv
    """
}
