/*
Calculate chromosome-specific stats
*/

params.options = [:]

process COMPUTE_DECONVOLUTION_STAT_CHR {
    tag "${sampleId}, ${chr}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"

    if (workflow.profile == 'aws') {
        label 'xlarge'
    } else {
        label 'cpu_medium'
        label 'memory_medium'
    }

  input:
    tuple val(sampleId), val(chr), path(bedpe), path(beads)
    val images_pulled

  output:
    tuple val(sampleId), val(chr), path("${sampleId}.*_overlapCount.csv.gz"), emit: overlap_count

  script:
    """
    atacFragOverlapMetricsChr.py \
      -nc 6 \
      -hq ${beads} \
      -be ${bedpe} \
      -csv ${sampleId}.${chr}_overlapCount.csv.gz \
      -c ${sampleId}.${chr}_ncCount.tsv \
      -rt 4 \
      -m ${params.options.mergeMethod.get(sampleId)} \
      -r ${params.options.rounding}
    """
}
