/*
Counting number of reads per gene annotation per cell
*/

params.options = [:]

process UMI_TOOLS_COUNT {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(bam), path(bam_index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.counts.tsv.gz'), emit: count

    script:
    include_introns = params.options.includeIntrons == true ? "^Unassigned|INTERGENIC" : "^Unassigned|INTRON|INTERGENIC"
    """
    mkdir -p tmp/

    umi_tools count \
        --method unique \
        --temp-dir "tmp/" \
        --per-cell \
        --per-gene \
        --extract-umi-method tag \
        --cell-tag XC \
        --umi-tag XM \
        --gene-tag XT \
        --skip-tags-regex "${include_introns}" \
        -I ${bam} \
        -S ${sampleId}.counts.tsv.gz
    """
}
