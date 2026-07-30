/*
Publish original and tidy-formatted pipeline parameters
*/

import org.yaml.snakeyaml.Yaml

process PUBLISH_PARAMETERS {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    if (!params.options.reportsDir) {
        reportsDir = "${params.core.reportsDir}"
    } else {
        reportsDir = "${params.options.reportsDir}"
    }
    if (workflow.profile =~ /(awsbatch|tower)/) {
        publishDir "${reportsDir}/", pattern:'*copy*', mode: 'copy', overwrite: true // groovylint-disable-line
    } else {
        publishDir "${reportsDir}/", pattern:'*copy*', mode: 'move', overwrite: true // groovylint-disable-line
    }
    label 'cpu_medium'
    label 'memory_xxsmall'

    input:
    val inputFile
    val images_pulled

    output:
    path "full_params.yaml", emit: params_yaml
    path "input_command.txt", emit: command_txt
    path "input_params.yaml", optional: true, emit: input_params_yaml
    path "full_params_copy.yaml", emit: params_yaml_copy
    path "input_command_copy.txt", emit: command_txt_copy
    path "input_params_copy.yaml", optional: true, emit: input_params_yaml_copy

    script:
    def abridgedParams = [:] << params

    List filterParams = ['preset', 'options', 'bead', 'cell', 'barcode',
    'quality-score', 'out-dir', 'results-dir', 'reports-dir', 'assets-dir',
    'max-cpus', 'max-memory', 'max-time', 'error-strategy', 'downsample']
    for (int i = 0; i < filterParams.size(); i++) {
        abridgedParams.remove(filterParams[i])
    }

    List omnitionGStringParams = ['outputDir', 'resultsDir', 'reportsDir',
    'assetsDir', 'modulesDir', 'workflowsDir']
    // Update all of the GString values to be strings under the omnition key
    for (int i = 0; i < omnitionGStringParams.size(); i++) {
        abridgedParams['omnition'][omnitionGStringParams[i]] = // groovylint-disable-line
        abridgedParams['omnition'][omnitionGStringParams[i]].toString() // groovylint-disable-line
    }

    // Update the GString parameters to be strings under the assay key
    List assayLevelGStringParams = ['resultsDir', // groovylint-disable-line
    'reportsDir', 'fastqFiles'] // groovylint-disable-line
    List referenceLevelGStringParams = ['fasta', 'gtf']
    List assays = ['rna', 'atac', 'catac']
    for (int i = 0; i < assays.size(); i++) {
        def assay = assays[i]
        // Make sure the assay is not null
        if (abridgedParams[assay] != null) {
            for (int j = 0; j < assayLevelGStringParams.size(); j++) {
                def assayParam = assayLevelGStringParams[j]
                abridgedParams[assay][assayParam] = abridgedParams[assay][assayParam].toString() // groovylint-disable-line
            }
        }
    }

    // Cannot nest under previous loop due to linter, too nested
    for (int i = 0; i < assays.size(); i++) { // groovylint-disable-line
        def assay = assays[i] // groovylint-disable-line
        // Make sure the assay is not null
        if (abridgedParams[assay] != null) {
            for (int j = 0; j < referenceLevelGStringParams.size(); j++) {
                def assayParam = referenceLevelGStringParams[j]
                if (abridgedParams[assay]['reference'][assayParam].getClass() == java.util.ArrayList) { // groovylint-disable-line
                    for (int k = 0; k < abridgedParams[assay]['reference'][assayParam].size(); k++) { // groovylint-disable-line
                        abridgedParams[assay]['reference'][assayParam][k] = abridgedParams[assay]['reference'][assayParam][k].toString() // groovylint-disable-line
                    } // groovylint-disable-line
                } else { // groovylint-disable-line
                    abridgedParams[assay]['reference'][assayParam] = abridgedParams[assay]['reference'][assayParam].toString() // groovylint-disable-line
                } // groovylint-disable-line
            } // groovylint-disable-line
        } // groovylint-disable-line
    } // groovylint-disable-line

    Yaml yamlBuilder = new Yaml()
    String formattedYamlParams = yamlBuilder.dump(abridgedParams)

    String formattedInputFile
    if (inputFile.getClass() == nextflow.config.ConfigMap) {
        formattedInputFile = yamlBuilder.dumpAsMap(inputFile)
    }

    if (!inputFile) {
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

    stub:
    /* groovylint-disable-next-line GStringExpressionWithinString */
    """
    # Create expected output files
    outputs=("full_params.yaml"
        "input_command.txt"
        "full_params_copy.yaml"
        "input_command_copy.txt"

    for output in "\${outputs[@]}"
    do
        touch \${output}
    done
    """
}
