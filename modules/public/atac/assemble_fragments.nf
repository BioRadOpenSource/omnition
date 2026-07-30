/*
Create chromosome-specific fragment files
*/

process ASSEMBLE_FRAGMENTS {
    tag "${sampleId}, ${chr}"
    container "${params.container_image.dbg_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), val(chr), path(bam), path(index), path(bead), path(blocklist)
    val images_pulled

    output:
    tuple val(sampleId), val(chr), path("${sampleId}.*.frag.bedpe.gz"), emit: assemble_fragments
    tuple val(sampleId), path("*_read_counts.csv"), emit: count
    tuple val(sampleId), path("${sampleId}.*.bead_counts.tsv"), emit: bead_counts

    script:
    """
    omnition atac assemble-fragments \
        --bam-file ${bam} \
        --blocklist ${blocklist} \
        --sample ${sampleId} \
        --contig ${chr} \
        --barcode-tag XB \
        --out ${sampleId}.${chr}.frag.bedpe.gz \
        --bead-counts ${sampleId}.${chr}.bead_counts.tsv \
        --read-counts ${sampleId}_${chr}_assemble_fragments_output_read_counts.csv \
        --mapping-quality ${params.options.mapQualityThreshold} \
        --max-insert-size ${params.options.maxInsertSize}
    """
}
