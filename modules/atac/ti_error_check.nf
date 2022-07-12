/*
Compile all TI errors
*/

params.options = [:]

process TI_ERROR_CHECK {
    if (workflow.profile == 'aws') {
        label 'small'
  } else {
        label 'cpu_xsmall'
        label 'memory_xxsmall'
    }
    errorStrategy 'terminate'
    publishDir "${params.reportsDir}", pattern: 'TI_run_errors.txt', mode: 'copy', overwrite: true

  input:
    path error_files
    val images_pulled

  output:
    path('TI_run_errors.txt')
    val true, emit:passed

  script:
    """
  > TI_run_errors.txt
  for FILE in ${error_files}; do
    cat \$FILE >> TI_run_errors.txt
  done
  error_check=\$(cat TI_run_errors.txt | wc -l)

  if [[ \$error_check -ne 0 ]]; then
    if ${params.options.tierroroverride}; then
      echo "Errors detected, but overridden"
    else
      echo "[ERROR] Error(s) detected for the following fastq files."
      cat TI_run_errors.txt
      exit 1
    fi
  fi
  """
}
