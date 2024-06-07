/*
Creating refFlat-formatted reference file
*/

params.options = [:]

process MAKE_REFFLAT {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path gtf
    val images_pulled

    output:
    path 'annotation.refflat', emit: refflat

    script:
    """
    # Creating extended format refFlat file
    gtfToGenePred -genePredExt ${gtf} annotation.refflat.tmp

    # Reformatting so the output file looks like the basic refFlat output from
    # gtfToGenePred but adds the gene name as the first column
    paste <(cut -f 12 annotation.refflat.tmp) <(cut -f 1-10 annotation.refflat.tmp) > annotation.refflat
    """
}
