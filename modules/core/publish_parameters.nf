/*
Publish original and tidy-formatted pipeline parameters
*/

import org.yaml.snakeyaml.Yaml

params.options = [:]

process PUBLISH_PARAMETERS {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile =~ /(awsbatch|tower)/) {
        publishDir "${params.options.reportsDir}", pattern:'*copy*', mode: 'copy', overwrite: true // groovylint-disable-line
    } else {
        publishDir "${params.options.reportsDir}", pattern:'*copy*', mode: 'move', overwrite: true // groovylint-disable-line
    }
    label 'cpu_medium'
    label 'memory_xxsmall'

    input:
    val inputFile
    val images_pulled

    output:
    path "full_params.yaml", emit: params_yaml
    path "input_command.txt", emit: commmand_txt
    path "input_params.yaml", optional: true, emit: input_params_yaml
    path "full_params_copy.yaml", emit: params_yaml_copy
    path "input_command_copy.txt", emit: commmand_txt_copy
    path "input_params_copy.yaml", optional: true, emit: input_params_yaml_copy

    script:
    def abridgedParams = [:] << params

    List filterParams = ['preset', 'options', 'bead', 'cell', 'barcode',
    'quality-score', 'out-dir', 'results-dir', 'reports-dir', 'assets-dir',
    'max-cpus', 'max-memory', 'max-time', 'error-strategy', 'downsample']
    for (int i = 0; i < filterParams.size(); i++) {
        abridgedParams.remove(filterParams[i])
    }

    Yaml yamlBuilder = new Yaml()
    String formattedYamlParams = yamlBuilder.dump(abridgedParams)

    // Adding condition for when params files are supplied via Seqera Platform
    String formattedInputFile
    if (inputFile.getClass() == nextflow.config.ConfigMap) {
        formattedInputFile = yamlBuilder.dumpAsMap(inputFile)
    }

    if (params.options.paramsFile == null) {
        """
        # Echo the pipeline command into a text file
        echo ${workflow.commandLine} > input_command.txt
        echo "${formattedYamlParams}" | yq | yq -Y >> full_params.yaml
        echo "No input parameters file supplied"

        cp full_params.yaml full_params_copy.yaml
        cp input_command.txt input_command_copy.txt
        """
    } else {
        """
        # Echo the pipeline command into a text file
        echo ${workflow.commandLine} > input_command.txt
        echo "${formattedYamlParams}" | yq | yq -Y >> full_params.yaml

        echo "${formattedInputFile}" >> input_params.yaml
        cp full_params.yaml full_params_copy.yaml
        cp input_command.txt input_command_copy.txt
        cp input_params.yaml input_params_copy.yaml
        """
    }
}
