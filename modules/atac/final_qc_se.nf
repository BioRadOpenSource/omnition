/*
Calculate final QC stats
*/

params.options = [:]

process FINAL_QC_SE {
    tag "${sampleId}"
    container "bioraddbg/omnition-r:${workflow.manifest.version}"

    if (workflow.profile == 'aws') {
        label 'xlarge'
  } else {
        label 'cpu_medium'
        label 'memory_xlarge'
    }

  input:
    tuple val(sampleId), path(fragments), path(index), path(qc), path(dict)
    path tss
    val images_pulled

  output:
    tuple val(sampleId), path("${sampleId}.QCstats.csv"), emit: qc_stats
    tuple val(sampleId), path("${sampleId}.QCstats.csv"), path(dict), emit: qc_compile

  script:
    """
  PEAK_FILE_GO="none"
  ONE_TO_ONE_GO="no"
  if [ "${params.options.mixed}" = "true" ]; then
    SPECIES_MIX_GO="yes"
  else
    SPECIES_MIX_GO="no"
  fi

  atacQcAdvancedSe.R \
    ${fragments} \
    ${tss} \
    \$PEAK_FILE_GO \
    ${qc} \
    \$SPECIES_MIX_GO \
    \$ONE_TO_ONE_GO \
    ${dict}
  """
}
