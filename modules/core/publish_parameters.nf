/*
Publish original and tidy-formatted pipeline parameters
*/

import org.yaml.snakeyaml.Yaml

process PUBLISH_PARAMETERS {
    container "bioraddbg/omnition-r:${workflow.manifest.version}"
    publishDir "${params.reportsDir}", mode: 'move', overwrite: true

    if (workflow.profile == 'aws') {
        label 'small'
    } else {
        label 'cpu_medium'
        label 'memory_xxsmall'
    }

    input:
    val images_pulled

    output:
    path "params.yaml", emit: params_yaml

    script:
    def abridgedParams = [:] << params
    abridgedParams.remove("preset")
    abridgedParams.remove("options")

    Yaml yamlBuilder = new Yaml()
    String formattedYamlParams = yamlBuilder.dump(abridgedParams)

    """
    echo "${formattedYamlParams}" >> params.yaml
    """
}
