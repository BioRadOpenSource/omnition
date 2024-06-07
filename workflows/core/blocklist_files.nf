workflow GET_BLOCKLIST_FILES {
    take:
    runParams

    main:
    if (workflow.profile =~ /(awsbatch|tower)/ && runParams.reference?.blocklist != null) {
        runParams.reference.blocklistFiles = [file("${runParams.reference.blocklist}")]
    }
}
