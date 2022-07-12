/*
Create chromosome-specific fragment files
*/

params.options = [:]

process ASSEMBLE_FRAGMENTS {
    tag "${sampleId}, ${chr}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"

    if (workflow.profile == 'aws') {
        label 'small'
    } else {
        label 'cpu_xsmall'
        label 'memory_xsmall'
    }

  input:
    tuple val(sampleId), val(chr), path(bam), path(index), path(bead)
    val images_pulled

  output:
    tuple val(sampleId), val(chr), path("${sampleId}.*.frag.bedpe.gz"), emit: assemble_fragments
    tuple val(sampleId), val(chr), path("*_read_counts.csv"), emit: stats

  script:
    """
    atacAssembleFragments.py \
        --input ${bam} \
        --sample ${sampleId} \
        --chromosome ${chr} \
        --mapping-quality ${params.options.qualityThreshold} \
        --insert-size ${params.options.maxInsertSize}
    """
}
