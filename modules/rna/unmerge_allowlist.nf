/*
Convert allowlist from drop barcodes to bead allowlist
*/

params.options = [:]

process UNMERGE_ALLOWLIST {
    tag "${sampleId}"
    container "bioraddbg/omnition-dbg:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(sampleId), path(allowlist), path(translate)
    val images_pulled

    output:
    tuple val(sampleId), path("*_unmerged_allowlist.csv"), emit: allowlist

    script:
    """
    rnaUnmergeAllowlist.py -a ${allowlist} -b ${translate}
    """
}
