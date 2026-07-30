/*
Split bam files by chromosome
*/

process SPLIT_BAM {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
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

    publicAtacNamesSplitFilt.py \
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
