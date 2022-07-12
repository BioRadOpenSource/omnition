/*
Check the read counts for each TI
*/

params.options = [:]

process CHECK_TI_COUNTS {
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
    path config
    val images_pulled

  output:
    path("*_error.txt"), emit: error_files
    path("*_reads.txt"), emit: split
    path("*_messages.txt"), emit: messages

  script:
    """
  touch "${sample_id}_error.txt"
  touch "${sample_id}_messages.txt"
  atacCountTiReads.py -i ${fastq[0]} -l $asset -c $config -s $sample_id -m ${sample_id}_messages.txt
  """
}
