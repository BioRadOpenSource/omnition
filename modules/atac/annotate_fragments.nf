/*
Annotate chromosome-specific fragment files
*/

params.options = [:]

process ANNOTATE_FRAGMENTS {
    tag "${sampleId}, ${chr}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), val(chr), path(bead), path(bedpe), path(blocklist)
    val images_pulled

    output:
    tuple val(sampleId), val(chr), path("${sampleId}.*.frag.bedpe.annotated.tsv.gz"), emit: annotate_fragments
    tuple val(sampleId), val(chr), path('*_read_counts.csv'), emit: stats
    tuple val(sampleId), path("${sampleId}.*.bead_counts.tsv"), emit: bead_counts

    script:
    """
    atacAnnotateFragments.R \
        ${blocklist} \
        ${bedpe} \
        ${bead} \
        ${sampleId}.${chr}.frag.bedpe.annotated.tsv \
        ${sampleId}.${chr}.bead_counts.tsv
    pigz -p ${task.cpus} ${sampleId}.${chr}.frag.bedpe.annotated.tsv

    READCOUNTFILE="${sampleId}_${chr}_annotate_fragments_output_read_counts.csv"
    printf "sample,process,metric,count\n" > \$READCOUNTFILE
    printf "${sampleId},annotate_fragments_${chr},input,\$(cat input_frags.tsv)\n" >> \$READCOUNTFILE
    printf "${sampleId},annotate_fragments_${chr},output,\$(cat output_frags.tsv)\n" >> \$READCOUNTFILE
    """
}
