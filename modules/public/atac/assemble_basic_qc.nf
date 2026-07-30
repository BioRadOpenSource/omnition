/*
Calculate basic QC stats
*/

process ASSEMBLE_BASIC_QC {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_xxsmall'

    input:
    tuple val(sampleId), path(stats)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.basicQC.tsv"), emit: assemble_basic_qc

    script:
    """
    cat \$(find . -name "*_frag.sumstats.tsv" -a -not -name \
        "${sampleId}.*${params.options.mitoContig}_frag.sumstats.tsv") \
        > ${sampleId}.basicQC-temp.tsv

    MITO_READS=\$(cat \$(find . -type f,l | \
        grep '${sampleId}.*${params.options.mitoContig}_frag.sumstats.tsv') | \
        wc -l)

    if [[ \$MITO_READS -gt 0 ]]; then
        cat \$(find . -name "${sampleId}.*${params.options.mitoContig}_frag.sumstats.tsv") > \
            ${sampleId}.basicQC-mito-temp.tsv
    else
        echo "NA\t0\t0\tNA" >> ${sampleId}.basicQC-mito-temp.tsv
    fi

    publicAtacQcSimple.R \
        ${sampleId}.basicQC-temp.tsv \
        ${sampleId}.basicQC-mito-temp.tsv \
        ${sampleId}.basicQC.tsv \
        ${params.options.mixed}
    """
}
