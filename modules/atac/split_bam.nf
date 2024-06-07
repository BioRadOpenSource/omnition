/*
Split bam files by chromosome
*/

params.options = [:]

process SPLIT_BAM {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_xlarge'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(bam), path(index), path(stats)
    path sizes
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.*.raw.bam"), path("${sampleId}.*.raw.bam.bai"),
        path("${sampleId}.*.read_bead.tsv.gz"), emit: split_bam optional true

    script:
    """
    MITO_CHR=\$(grep "${params.options.mitoContig}" ${sizes} | cut -f1 | tr "\\n" ,| awk 'BEGIN {FS=OFS=","} NF--')

    atacNamesSplitFilt.py \
        --input ${bam} \
        --name ${sampleId} \
        --output ./ \
        --barcode-tag XB \
        --bedtools-reference-genome ${sizes} \
        --mito-chr \$MITO_CHR \
        --ncores ${task.cpus} \
        --mapq ${params.options.mapQualityThreshold} \
        --max-insert ${params.options.maxInsertSize}
    """
}
