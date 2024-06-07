/*
Decompressing files
*/

params.options = [:]

process GUNZIP {
    tag "${file}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_xxsmall'

    input:
    path file
    val images_pulled

    output:
    path "${out}", emit: decompressed

    script:
    if (file.name =~ /gz$/) {
        out = file.baseName
        """
        gunzip -f ${file}
        """
    } else {
        out = file
        """
        touch ${file}
        """
    }
}
