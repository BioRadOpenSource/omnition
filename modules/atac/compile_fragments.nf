/*
Compile fragments by sample
*/
process COMPILE_FRAGMENTS {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.resultsDir}/${sampleId}/fragments", mode: "copy", overwrite: true

    if (workflow.profile == 'aws') {
        label 'large'
    } else {
        label 'cpu_large'
        label 'memory_medium'
    }

    input:
        tuple val(sampleId), path(fragments), path(index)
        val images_pulled

    output:
        tuple val(sampleId), path("${sampleId}.fragments.tsv.gz"),
        path("${sampleId}.fragments.tsv.gz.tbi"), emit:fragments

    script:
        """
        # Using `ls` to sort file names for consistency
        zcat \$(ls ${fragments}) | sort -V -k1,1 -k2,2 | bgzip > ${sampleId}.fragments.tsv.gz
        tabix -p bed ${sampleId}.fragments.tsv.gz
        """
}
