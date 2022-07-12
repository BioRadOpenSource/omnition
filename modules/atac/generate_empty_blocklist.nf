/*
Create an empty blocklist
*/

params.options = [:]

process GENERATE_EMPTY_BLOCKLIST {
    if (workflow.profile == 'aws') {
        label 'small'
  } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }

  input:
    val images_pulled

  output:
    path 'blocklist.bed', emit:bed

  script:
    '''
  touch blocklist.bed
  '''
}
