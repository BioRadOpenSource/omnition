/*
Reannotate fragments files based on merging scheme
*/

params.options = [:]

process REANNOTATE_FRAGMENTS {
    tag "${sampleId}, ${chr}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"

    if (workflow.profile == 'aws') {
        label 'xlarge'
  } else {
        label 'cpu_medium'
        label 'memory_medium'
    }

  input:
    tuple val(sampleId), val(chr), path(bedpe), path(dict)
    val images_pulled

  output:
    tuple val(sampleId), val(chr), path("${sampleId}.*.frag.bedpe.annotated.dedup.tsv"), emit: chr_reannotate_fragments
    tuple val(sampleId), path("${sampleId}.*.frag.bedpe.annotated.dedup.tsv"), emit: reannotate_fragments
    tuple val(sampleId), path("${sampleId}.*_frag.sumstats.tsv"), emit: frag_sumstats
    tuple val(sampleId), val(chr), path('*_read_counts.csv'), emit: stats

  script:
    """
  atacReannotateFragments.R \
    ${bedpe} \
    ${dict} \
    ${sampleId}.${chr}.frag.bedpe.annotated.dedup.tsv \
    ${sampleId}.${chr}_frag.sumstats.tsv

  READCOUNTFILE="${sampleId}_${chr}_reannotate_fragments_output_read_counts.csv"
  printf "sample,process,metric,count\n" > \$READCOUNTFILE
  printf "${sampleId},reannotate_fragments_${chr},input,\$(cat input_frags.tsv)\n" >> \$READCOUNTFILE
  printf "${sampleId},reannotate_fragments_${chr},output,\$(cat output_frags.tsv)\n" >> \$READCOUNTFILE
  """
}
