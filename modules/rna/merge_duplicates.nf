/*
Convert duplicate count file from bead barcodes to drop barcodes
*/

params.options = [:]

process MERGE_DUPLICATES {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_large'

    input:
    tuple val(sampleId), path(duplicates), path(translate)
    val images_pulled

    output:
    tuple val(sampleId), path("*.merged_duplicate_counts.csv"), emit: duplicate_count

    script:
    """
    rnaMergeDuplicates.py -d ${duplicates} -b ${translate}
    """
}
