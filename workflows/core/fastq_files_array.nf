workflow GET_FASTQ_FILES {
    take:
    runParams

    main:
    if (workflow.profile =~ /(awsbatch|tower)/) {
        runParams.fastqFiles = file("${runParams.input}" + Core.fastqRegEx()).name.toArray()
    } else {
        runParams.fastqFiles = new File(runParams.input).listFiles()
    }
}
