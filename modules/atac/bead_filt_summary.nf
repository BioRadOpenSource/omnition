/*
Bead filtration summary table
*/

params.options = [:]

process BEAD_FILT_SUMMARY {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'large'
  } else {
        label 'cpu_small'
        label 'memory_xxsmall'
    }

  input:
    tuple val(sampleId), path(stats), path(quant), path(params)
    path counts, stageAs: 'fastqTIreadcounts.csv'
    val images_pulled

  output:
    tuple val(sampleId), path("${sampleId}.beadFiltSummary.csv"), emit:index_summary
    tuple val(sampleId), path("${sampleId}.sampleBeadFiltSummary.csv"), emit:sample_summary

  script:
    """
  atacBeadFiltSummary.py -s ${sampleId} -i ./
  """
}
