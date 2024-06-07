/*
Compute transcription start site matrix
*/

params.options = [:]

process COMPUTE_TSS_MATRIX {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/bulkQC", pattern: '*.tss_data_matrix.gz', mode: 'copy',
        overwrite: true
    label 'cpu_small'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bam), path(bai), path(reference_tss)
    val images_pulled

    output:
    tuple val(sampleId), path('*.tss_data_matrix.gz'), emit: tss_matrix // groovylint-disable-line

    script:
    """
    # Computing the TSS window as the entire range, plus the central TSS base
    FULL_WINDOW=\$(( ${params.options.tssWindowSize} + 1 ))

    # Get the order of contigs in the alignments
    samtools view -H ${bam} | grep SQ | cut -f 2,3 | awk '{ sub(/^SN:/, ""); sub(/LN:/, ""); print;}' > genome.txt

    # Sort reference TSS to match BAM order and shift window by one base
    cut -f1 genome.txt |
        while read -r line;
        do
            grep -w "\$line" ${reference_tss} |
            awk -v WID=\$(( \$FULL_WINDOW - 1 )) 'BEGIN {OFS="\t"} {
                if(\$3 - \$2 == WID )
                {
                print \$1,\$2-1,\$3,\$4
                }
            }' >> tmp.tss.bed ;
        done

    # Calculate coverage and format into a matrix
    bedtools coverage -a tmp.tss.bed -b ${bam} -d -sorted -g genome.txt | \
        cut -f6 | \
        awk -v n=\$FULL_WINDOW '{ row = row \$1 "\\t"; if (NR % n == 0) { print row; row = "" } }' | \
        cut -f1-\$FULL_WINDOW |
        pigz -p ${task.cpus} > ${sampleId}.tss_data_matrix.gz

    """
}
