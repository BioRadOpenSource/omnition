/*
Caluclate the fraction of reads in peaks
*/

params.options = [:]

process FRACTION_OF_READS_IN_PEAKS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bam), path(index), path(bed)
    val images_pulled

    output:
    tuple val(sampleId),  path('*.frip.csv'), emit: frip

    script:
    """
    bedtools intersect -a ${bam} -b ${bed} > readsInPeaks.bam
    readsInPeaks=\$(samtools view -c readsInPeaks.bam)
    totalReads=\$(samtools view -c ${bam})
    frip=\$(printf "%.2f\n" \$((10**2 * \$readsInPeaks/\$totalReads))e-2)

    printf "sample,process,metric,value\n" > ${sampleId}.frip.csv
    printf "${sampleId},frip,readsInPeaks,\${readsInPeaks}\n" >> ${sampleId}.frip.csv
    printf "${sampleId},frip,totalDedupedReads,\${totalReads}\n" >> ${sampleId}.frip.csv
    printf "${sampleId},frip,fractionOfReadsInPeaks,\${frip}\n" >> ${sampleId}.frip.csv
    """
}
