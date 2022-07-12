/*
Create a json config file of TIs for DEAD
*/

params.options = [:]

process TI_DEAD_CONFIG {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'small'
  } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }
    errorStrategy 'terminate'

  input:
    path asset
    path json
    val images_pulled

  output:
    path("TI.json"), emit: config
    path("TIlen.txt"), emit: tilen

  script:
    """
  atacConfigToDeadJson.py -i ${asset} -m 1 -r ${params.options.tiread} -t ${params.options.i7asti}
  """
}
