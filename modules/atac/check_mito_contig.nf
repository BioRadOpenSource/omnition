/*
Check the mitochondrial contig exists in the reference files
*/

params.options = [:]

process CHECK_MITO_CONTIG {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path sizes
    val images_pulled

    output:

    script:
    """
    MITO_CHR=\$(grep "${params.options.mitoContig}" ${sizes} | cut -f1 | tr "\\n" ,| awk 'BEGIN {FS=OFS=","} NF--')
    # Ensure that the mitochondrial contig given isn't missing from the reference
    if [ -z "\$MITO_CHR" ] || [ "\$MITO_CHR" == " " ]; then
        echo \
        "[ERROR] Provided contig ${params.options.mitoContig} not found in reference files.  Please check reference. \
            Aborting pipeline" && exit 1
    fi

    """
}
