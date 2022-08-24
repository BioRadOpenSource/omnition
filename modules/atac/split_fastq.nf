/*
Split fastqs by TI
*/

params.options = [:]

process SPLIT_FASTQ {
    tag "${sample_id}-${index}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'small'
  } else {
        label 'cpu_medium'
        label 'memory_small'
    }
    errorStrategy 'terminate'

  input:
    tuple val(sample_id), val(index), path(fastq)
    val ti_checks_passed
    val images_pulled

  output:
    tuple val("${sample_id}-${index}"), path("${sample_id}-*.complete_debarcoded.split.fastq.gz"), emit: split

  script:
    """
  atacSplitFastq.py -r1 ${fastq[0]} -r2 ${fastq[1]} -i $index
  pigz -p ${task.cpus} ${sample_id}-${index}_R1.complete_debarcoded.split.fastq
  pigz -p ${task.cpus} ${sample_id}-${index}_R2.complete_debarcoded.split.fastq
  """
}
