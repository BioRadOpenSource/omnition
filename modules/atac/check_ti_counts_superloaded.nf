/*
Check the read counts for each TI in a superloaded sample
*/

params.options = [:]

process CHECK_TI_COUNTS_SUPERLOADED {
    tag "${sample_id}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'small'
  } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }
    errorStrategy 'terminate'

  input:
    tuple val(sample_id), path(fastq)
    path asset
    val images_pulled

  output:
    path("*_reads.txt"), emit: split
    path("*_messages.txt"), emit: messages
    val true, emit:passed

  script:
    """
  touch "${sample_id}_messages.txt"
  atacCountTiReadsSuperloaded.py -i ${fastq[0]} -l $asset -s $sample_id -m ${sample_id}_messages.txt
  """
}
