/*
Handle creation of TSS windows
*/

params.options = [:]

process GENERATE_TSS_WINDOWS {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}", mode: 'copy', overwrite: true
    if (workflow.profile == 'aws') {
        label 'small'
  } else {
        label 'cpu_medium'
        label 'memory_xxsmall'
    }

  input:
    path gtf
    path fasta
    path sizes
    val images_pulled

  output:
    path "TSS.${params.options.tssWindowSize}.bed", emit:tss

  script:
    """
  # Computing the TSS window as a distance from the central TSS base
  HALF_WINDOW=\$(( ${params.options.tssWindowSize} / 2 ))

  # Check for ENCODE basic tag
  BASIC_TAG="\$( grep -m 1 "tag \"basic\"" ${gtf} || true )"

  # Convert TSS to BED windows
  if [[ -z \$BASIC_TAG ]]; then
    awk 'BEGIN{OFS="\\t"} {if(\$3 == "transcript") {if(\$7 == "-") {v = \$5} \
      else {v = \$4} print \$1,v,v,\$7}}' ${gtf} | sortBed -i stdin | uniq > temp.TSS.bed
    bedtools slop -i temp.TSS.bed -g ${sizes} -b \${HALF_WINDOW} > TSS.${params.options.tssWindowSize}.bed
  else
    grep ${gtf} "tag \"basic\"" | awk 'BEGIN{OFS="\\t"} {if(\$3 == "transcript") {if(\$7 == "-") {v = \$5} \
      else {v = \$4} print \$1,v,v,\$7}}' | sortBed -i stdin | uniq > temp.TSS.bed
    bedtools slop -i temp.TSS.bed -g ${sizes} -b \${HALF_WINDOW} > TSS.${params.options.tssWindowSize}.bed
  fi
  """
}
