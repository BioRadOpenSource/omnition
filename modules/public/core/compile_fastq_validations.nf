/*
Sort error files into passed and failed
Halt the pipeline if any samples failed
*/

process COMPILE_FASTQ_VALIDATIONS {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}/", mode: 'copy', overwrite: true
    errorStrategy 'terminate'
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path files
    val images_pulled

    output:
    path("*_errors.txt"), includeInputs: true, emit: val_pass

    script:
    List passed = []
    List failed = []
    for (int i = 0; i < files.size(); i++) {
        if (files[i].size() == 0) {
            passed.add(files[i])
        }
        else {
            failed.add(files[i])
        }
    }

    failed.sort { it.name }

    // Create an output file to write the contents of the failed samples to
    def output = "fastq_validation_errors.txt"
    // Write the contents of the files in the failed list to the output file
    if (failed.size() != 0) {
        def outFile = new File("${params.core.reportsDir}${output}")
        outFile.withWriter { writer ->
            for (file in failed) {
                file.eachLine { line ->
                    writer.writeLine(line)
                }
            }
        }
        log.error("[$params.core.assay] FASTQ read files failed validation. " \
                + "Please see ${params.core.reportsDir}fastq_validation_errors.txt " \
                + "for more details.")
        """
        exit 1
        """
    } else {
        """
        # Touch the files to make them valid outputs
        for file in ${files}
        do
            touch \$file
        done
        """
    }
}
