/*
Count features per barcode
*/

params.options = [:]

process BEAD_COUNTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_small'
    label 'memory_large'

    input:
    tuple val(sampleId), path(bam), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.sequence_counts.csv'), emit: sequence_counts

    script:
    """
    # Generate a per-bead per-contig count table
    rnaUmiCounter.py -b ${bam} -bt XC -ut XM -c ${task.cpus}
    """
}
