/*
Calculating crosstalk statistics, etc. for mixed species experiments
*/

process SUMMARIZE_MIXED_SPECIES {
    tag "${sampleId}"
    container "${params.container_image.r_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/countMatrix/",
        pattern:'*.species_mix_counts.{rds,csv}', mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    tuple val(species1), val(species2)
    tuple val(sampleId), path(counts), path(allowlist)
    val images_pulled

    output:
    tuple val(sampleId), path("${sampleId}.species_mix_counts.rds"), emit: count
    tuple val(sampleId), path("${sampleId}.species_mix_counts.csv"), emit: count_csv
    tuple val(species1), val(sampleId), path("${sampleId}.${species1}.allowlist.csv"), emit: s1_allowlist
    tuple val(species2), val(sampleId), path("${sampleId}.${species2}.allowlist.csv"), emit: s2_allowlist
    tuple val(sampleId), path("${sampleId}_crosstalk_density.csv"), emit: crosstalk_density
    tuple val(sampleId), path("${sampleId}.crosstalk.csv"), emit: metrics

    script:
    """
    publicCoreSpeciesMixStats.R -p ${sampleId} -a ${allowlist} -ct ${params.options.crosstalkThreshold} \
        ${counts} ${species1} ${species2}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=8

    # Check that we have expected number of sample input files
    check_input_files "SUMMARIZE_MIXED_SPECIES" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}.species_mix_counts.rds" "${sampleId}.${species1}.allowlist.csv"
    "${sampleId}.${species2}.allowlist.csv" "${sampleId}_crosstalk_density.csv"
    "${sampleId}.species_mix_counts.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("sampleId" "species1" "allowlist" "species2" "workflow.manifest.version"
    "params.options.crosstalkThreshold" "counts" "params.options.shadowThreshold"
    "params.options.resultsDir")

    parameters=("${sampleId}" "${species1}" "${allowlist}" "${species2}"
    "${workflow.manifest.version}" "${params.options.crosstalkThreshold}" "${counts}"
    "${params.options.shadowThreshold}" "${params.options.resultsDir}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "SUMMARIZE_MIXED_SPECIES" "\$param_name" "\$parameter"
    done
    """
}
