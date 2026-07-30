/*
Aggregating all sample metrics together into a tidy format
*/

process AGGREGATE_METRICS {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}/", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_xsmall'
    if (!(workflow.profile =~ /(awsbatch|tower)/)) {
        beforeScript 'ulimit -Ss unlimited'
    }

    input:
    path input
    val assay
    val images_pulled

    output:
    path 'metric_summary*.csv', emit: summary

    script:
    if (input.name == "paths.txt") {
        """
        cat ${input} | parallel -j ${task.cpus} aws s3 cp s3:/{} ./
        python -m bragg -i ./ -w agg_met_work -o ./ -a "${assay}"
        """
    } else {
        """
        python -m bragg -i ./ -w agg_met_work -o ./ -a "${assay}"
        """
    }

    stub:
    """
    # Create expected output files
    outputs=("metric_summary.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.reportsDir" "workflow.manifest.version"
        "assay")

    parameters=("${params.options.reportsDir}" "${workflow.manifest.version}"
        "${assay}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "AGGREGATE_METRICS" "\$param_name" "\$parameter"
    done
    """
}

