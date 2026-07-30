
class PublicCore {

    /* groovylint-disable */
    // Function to set ANSI color codes
    static Map ansiColors(Boolean monochrome_logs) {
        Map codes = [:]

        codes['reset']  = monochrome_logs ? '' : "\033[0m"
        codes['dim']    = monochrome_logs ? '' : "\033[2m"
        codes['black']  = monochrome_logs ? '' : "\033[0;30m"
        codes['red']    = monochrome_logs ? '' : "\033[0;31m"
        codes['green']  = monochrome_logs ? '' : "\033[0;32m"
        codes['yellow'] = monochrome_logs ? '' : "\033[0;33m"
        codes['blue']   = monochrome_logs ? '' : "\033[0;34m"
        codes['purple'] = monochrome_logs ? '' : "\033[0;35m"
        codes['cyan']   = monochrome_logs ? '' : "\033[0;36m"
        codes['white']  = monochrome_logs ? '' : "\033[0;37m"

        return codes
    }

    // Function to create pipeline logo
    static String bioradLogo(monochrome_logs) {
        Map colors = ansiColors(monochrome_logs)
        def logo = ''

        logo = String.format(
            """\n
            -${colors.dim}-------------------------------------------------------------------------------------${colors.reset}-

                      ${colors.green}###################################################################${colors.reset}
                   ${colors.green}#########################################################################${colors.reset}
                 ${colors.green}#############################################################################${colors.reset}
                ${colors.green}#########          #    ##         ######          ####     ###         #######${colors.reset}
                ${colors.green}#########   ###   ##   #    ####    #####   ####   ###      ###   ###    ######${colors.reset}
                ${colors.green}########          ##   #   #####    #             ###   ##   #   #####   ######${colors.reset}
                ${colors.green}########   ####   #    #   ####    #####   ####   #          #   ####   #######${colors.reset}
                ${colors.green}#######          ##   ###        #######   ####      #####             ########${colors.reset}
                 ${colors.green}#############################################################################${colors.reset}
                   ${colors.green}#########################################################################${colors.reset}
                     ${colors.green}###################################################################${colors.reset}
                                           ${colors.white}Omnition Analysis Software${colors.reset}

            -${colors.dim}-------------------------------------------------------------------------------------${colors.reset}-
            \n""".stripIndent()
       )

        return logo
    }
    /* groovylint-enable */

    static String fastqRegEx() {
        return "*_R{1,2}_*.{fq,fastq,fq.gz,fastq.gz}"
    }

    // Function to create help message
    static String helpMessage(monochrome_logs) { // groovylint-disable-line
        // logo = bioradLogo(monochrome_logs)

        def message = '' // groovylint-disable-line

        message = String.format(
            '''\n
            Standard usage
            nextflow run BioRadOpenSource/omnition -params-file <path to parameters file>

            Nextflow arguments:
            -params-file            Path to the parameter definition file
            -profile                Execute using Docker containers or Singularity containers(default: Singularity).
            -resume                 Resume a previous pipeline execution using the cached results(default: no resume).
            -r, -revision           Specify git tag, branch, or commit hash from which

            Bio-Rad arguments:
            --errorStrategy         Pipeline behavior when a process fails(default: ignore).
            --force                 Overwrite data in output directory(default: no overwrite).
            --help                  Show this message
            \n'''.stripIndent()
       )

        return message
    }

    static String nameFormatMessage() {
        return "FASTQ read files are present in the input directory, " \
                    + "but the name format is incorrect. Check parameters and/or see " \
                    + "documentation for file naming guidelines."
    }

    static String fileNameRegEx() {
        return '[a-zA-Z0-9_]+_S[0-9]+(_L[0-9]{3})?_R[0-9]+_001.f(ast)?q(.gz)?'
    }

    // Function to log parameter errors
    static String paramError(assay, key, log) {
        log.error("ERROR: [$assay] Must set '$key' parameter.")
        System.exit(1)
        return
    }

    // Function for counting the number of elements in an object
    static int countElement(object) {
        int value

        if (object.class == String || object.class == org.codehaus.groovy.runtime.GStringImpl) {
            value = 1
        } else if (object.class == ArrayList) {
            value = object.size()
        }

        return value
    }

    // Function to check 'output' parameter
    static String validateOutput(params, coreParams, log) {
        String value
        File file

        // Check if a specific output directory is set
        if (coreParams?.outputDir != null) {
            value = coreParams.outputDir
            if (value.endsWith("/")) {
                value = value.substring(0, value.length() - 1) // remove trailing '/'
            }
            file = new File(value)
        } else {
            value = params.omnition.outputDir // this is the default "./results" dir
            file = new File(value)
        }
        if (params.force == null) {
            // Check if the specified directory exists and if has anything in it other than the trace file
            if ((file.list() != null) && (file.list().length > 0) && !onlyTraceFileInOutDir(file, params)) {
                // Check if the force flag was not supplied
                log.error("ERROR: The specified output directory is not empty. "
                    + "This can be overridden with the --force flag: $value.")
                System.exit(1)
            } else {
                return value
            }
        } else {
            return value
        }
        return null
    }

    // Function to validate the path to the parameter file
    static <Type> Type validateParamsFile(params) {
        Type value

        // Check if params.omnition.paramsFile is populated (not null and not empty)
        if (params.omnition?.paramsFile) {
            value = params.omnition.paramsFile
        } else {
            value = null
        }
        return value
    }

    // Function to check if prefix is set, and if it is, override default
    static String validatePrefix(params, assayParams, log) {
        String value

        // Check if set
        if (params.prefix != null && params.prefix.class != String) {
            log.error("ERROR: [$assayParams.assay] Prefix parameter must be a string.")
            System.exit(1)
            return
        } else if (params.prefix != null) {
            value = params.prefix.concat('-')
            return value
        }
        return "" // groovylint-disable-line
    }

    // Function to check the reference 'directory' parameter
    static String validateReferenceDirectory(params, log) {
        String value

        // Check if set
        if (params.reference?.directory != null) {
            value = params.reference.directory
            return value
        }
        log.error("ERROR: [$params.assay] Must set the reference 'directory' parameter.")
        System.exit(1)
        return
    }

    // Function to check the 'workflow' parameter
    static String validateWorkflow(params, log) {
        String value

        // Check if set
        if (params.workflow != null) {
            value = params.workflow
            // Check if provided value is valid
            if (!(value in [ 'reference', 'analysis', 'full' ])) {
                log.error("ERROR: [$params.assay] The 'workflow' parameter must be \
                    'reference', 'analysis', or 'full'.")
                System.exit(1)
            }
            return value
        }
        PublicCore.paramError(params.assay, 'workflow', log)
        return
    }

    // Function to ensure the input file names match the Illumina naming convention
    // groovylint-disable
    // https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm 
    // groovylint-enable
    static boolean validateInputNames(directory) {
        def filename = fileNameRegEx()
        def validInputs = directory.listFiles().findAll {
                    it.name.matches(filename) }
        // checks for .fastq.gz or .fq.gz mispelling
        // skips directories
        // skips hidden files and files without extentions
        def files = directory.listFiles().stream().findAll {
                    !it.directory &&
                    it.name.indexOf('.') > 0 &&
                    it.name.matches(it.name.take(it.name.indexOf('.')) + "*.[a-z]{1,6}(.gz)?") }
        return files.size() == validInputs.size()
    }

    static boolean validateInputNamesAWS(sampleFiles) {
        def filename = fileNameRegEx()

        int awsFileListCount = sampleFiles
            .stream()
            .findAll { it.matches(filename) }
            .collect()
            .size()

        int filesCount = sampleFiles.stream().count()
        return filesCount == awsFileListCount
    }

    static String validateInputAWS(params, sampleFiles, log, messages) {
        if (sampleFiles != null) {
            boolean test = sampleFiles.stream().anyMatch {
                    it.toString().matches('.*\\/Undetermined[\\w\\.\\-]+') } // groovylint-disable-line
            // Warn if FASTQ files beginning with 'Undetermined' are found
            if (test) {
                log.warn("[$params.assay] Undetermined FASTQ read files will be ignored.")
                messages.add("WARN: [$params.assay] Undetermined FASTQ read files will be ignored.")
            }
            if (sampleFiles.size() != 0) {
                if (!PublicCore.validateInputNamesAWS(sampleFiles)) {
                    log.error("ERROR: [$params.assay] All FASTQ files in the 'input' parameter S3 bucket must "\
                      + "match the Illumina naming convention SampleName_S1_L001_R1_001.fastq.gz " \
                      + "OR SampleName_S1_R1_001.fastq.gz")
                    System.exit(1)
                }
                return params.input
            }
            log.error("ERROR: [$params.assay] The 'input' parameter S3 bucket does not contain any files")
            System.exit(1)
        } else {
            PublicCore.paramError(params.assay, 'input', log) // groovylint-disable-line
        }
    }

    // Function to check 'input' directory parameter
    static String validateInput(params, log, messages) {
        String value
        File file

        // Check if set
        if (params.input != null) {
            value = params.input
            file = new File(value)
            // Check if directory exists
            if (file.exists()) {
                boolean test = file.listFiles().stream().anyMatch {
                    it.toString().matches('.*\\/Undetermined[\\w\\.\\-]+') } // groovylint-disable-line
                // Warn if FASTQ files beginning with 'Undetermined' are found
                if (test) {
                    log.warn("[$params.assay] Undetermined FASTQ read files will be ignored.")
                    messages.add("WARN: [$params.assay] Undetermined FASTQ read files will be ignored.")
                }
                /* groovylint-disable */
                if (!PublicCore.validateInputNames(file)) {
                    log.error("ERROR: [$params.assay] All FASTQ files in the 'input' parameter directory must "\
                      + "match the Illumina naming convention SampleName_S1_L001_R1_001.fastq.gz " \
                      + "OR SampleName_S1_R1_001.fastq.gz") 
                    System.exit(1)
                }
                return value
            }
            log.error("ERROR: [$params.assay] The 'input' parameter directory does not exist: $value") 
            System.exit(1)
        } else {
            PublicCore.paramError(params.assay, 'input', log)             
        }
        /* groovylint-enable */
    }

    // Function to ensure that species name / prefix match between fasta-gtf pairs
    static void validateReferenceNames(assayParams, log) {
        // Get just the prefix from the file names
        List gtfSpeciesList = assayParams.reference.gtf
            .stream()
            .filter(Objects::nonNull)
            .map(i -> i - ~/^\/.*\//)
            .map(i -> i - ~/\..*$/)
            .distinct()
            .collect()

        List fastaSpeciesList = assayParams.reference.fasta
            .stream()
            .filter(Objects::nonNull)
            .map(i -> i - ~/^\/.*\//)
            .map(i -> i - ~/\..*$/)
            .distinct()
            .collect()

        // The prefixes should match and be in the same order
        if (gtfSpeciesList == (fastaSpeciesList)) {
            log.info("INFO: [$assayParams.assay] FASTA and GTF files match.")
        } else {
            log.error("ERROR: [$assayParams.assay] FASTA and GTF file prefixes do not match.")
            System.exit(1)
        }
    }

    static List<String> getSpeciesSchema(assay, gtf, log) {
        // Need to split on comma or space separation
        List speciesList = gtf.split(/[\s,]+/)
        .stream()
        .filter(Objects::nonNull)
        .map(i -> i - ~/^s3:\/\/.*\//)
        .map(i -> i - ~/^\/.*\//)
        .map(i -> i - ~/\..*$/)
        .distinct()
        .collect()

        // Checking if files were found
        if (speciesList.empty) {
            log.error("ERROR: [$assay] No GTF files found in input directory. Check"
                    + "parameters and file name requirements.") // groovylint-disable-line
            System.exit(1)
        } else {
            return speciesList
        }
    }

    // Function for extracting all sample IDs matching a given pattern (NOTE: filters out "Undetermined" samples)
    static List<String> getSampleIds(assay, fileList, pattern, log, messages) {
        // Formatting inputs
        List globList = pattern instanceof List ? pattern : [pattern]

        // Creating list of all file prefixes
        List prefixList = []
        for (int i = 0; i < globList.size(); i++) {
            List matchList = fileList
                .stream()
                .map(j -> getPrefix(j.toString(), globList[i].toString(), log))
                .collect()

            prefixList.addAll(matchList)
        }

        // Creating list of sample IDs for all files matching all globs
        List sampleIdList = prefixList
            .stream()
            .filter(Objects::nonNull)
            .filter(i -> !i.startsWith('Undetermined'))
            .map(i -> i - ~/(_L[0-9][0-9][0-9])$/)
            .distinct()
            .collect()

        // Checking if files were found
        if (sampleIdList.empty) {
            log.error("ERROR: [$assay] No FASTQ files found in input directory. "
                + "Check parameters and file name requirements.") // groovylint-disable-line
            System.exit(1)
        } else {
            log.info("INFO: [$assay] Sample IDs: $sampleIdList")
            messages.add("INFO: [$assay] Sample IDs: $sampleIdList")
            return sampleIdList
        }
    }

    // Function for parsing sample bead and cell settings into a CSV format
    static List<String> getSampleSettings(params) {
        // Check if this is a simple settings map (format, length) or a sample-keyed map
        def firstEntry = params.entrySet().iterator().next()
        println("First entry in params: ${firstEntry.key} = ${firstEntry.value} (${firstEntry.value?.getClass()})")

        if (firstEntry.value instanceof String || firstEntry.value instanceof Number) {
            println("WARNING: getSampleSettings received simple key-value pairs, not sample-keyed settings")
            println("Expected structure: [sampleId: [setting: value]], got: [setting: value]")
            return []
        }

        List sampleParamsList = params
            .entrySet()
            .stream()
            .map { s1 ->
                String setting = s1.key
                List result = s1.value
                    .entrySet()
                    .stream()
                    .map { s2 ->
                        String sample = s2.key
                        String value = s2.value
                        return [ sample, setting, value ]
                    }
                    .collect()
                return result
            }
            .collect()

        return sampleParamsList
    }

    // Function to report a warning if the number of samples exceeds 12
    static <Type> Type validateNumberOfSamples(assayParams, messages) {
        // Check if set
        if (assayParams.sampleIds != null) {
            int value = assayParams.sampleIds.size()
            // Check if class type is correct
            if (value > 12) {
                messages.add("INFO: [$assayParams.assay] Over 12 samples are being analyzed, "
                + "this exceeds the design limit for the report. Output files will be unaffected.")
            }
        }
    }

    // Function to ensure the user-provided merged sample names match the Illumina naming convention
    // groovylint-disable
    //https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm 
    // groovylint-enable
    static boolean validateConcatenateSampleName(value) {
        def sample = '[a-zA-Z0-9_]+_S[0-9]+'
        return value ==~ sample
    }

    // Funtion to validate which files to concatenate and what to label the concatenated sample
    static <Type> Type updateSampleIds(assayParams) {
        // Initializing variables
        List sampleIdList = []
        sampleIdList = assayParams.sampleIds
        return sampleIdList
    }

    // Function for checking if the only file in the output directory is a trace file
    private static boolean onlyTraceFileInOutDir(dir, params) {
        // Gather information on the output dir files and count differences between
        // expected files
        List expected = [".nextflow", ".nextflow.log", "pipeline_info", "work"]
        int diff = (dir.list() - expected).size()

        if (diff > 0 || dir.list().findIndexValues { it ==~ "pipeline_info" } == null) { // groovylint-disable-line
            return false
        }

        // Get the pipeline info subDir to check for only trace file existing
        int index = dir.list().findIndexValues { it ==~ "pipeline_info" } [0] // groovylint-disable-line
        File subDir = dir.listFiles()[index]

        // Get trace name string that allows for setting prefix
        String trace = params.prefix ? (params.prefix + "-omnition-trace-") : ("omnition-trace-")

        if (!subDir.directory) {
            return false
        }
        if (subDir.list().length > 1) {
            return false
        }
        return subDir.listFiles()[0].toString().contains(trace)
    }

    // shared function to load core preset params
    private static <Type> void loadCorePresetParams(params) {
        def presetsCore = [:]
        for (key in params.preset.core.keySet()) {
            if ((key in params.core.keySet()) && (params.preset.core[key] instanceof Map)) {
                params.preset.core[key].keySet().findAll {
                    !(it in params.core[key].keySet())
                }.each {
                    presetsCore[key] = presetsCore.get(key, [:]); presetsCore[key][it] = params.preset.core[key][it]
                }
            } else {
                if (!(key in params.core.keySet())) {
                    presetsCore[key] = params.preset.core[key]
                }
            }
        }

        // add the presetsCore to the params.core, but if the preset is a map, merge it with the existing map
        for (key in presetsCore.keySet()) {
            if (key in params.core.keySet() && params.core[key] instanceof Map) {
                params.core[key].putAll(presetsCore[key])
            } else {
                params.core[key] = presetsCore[key]
            }
        }
        return
    }

    // shared function to load do_merging preset params
    private static <Type> void loadDOMergingPresetParams(params) {
        def presetsDoMerging = [:]
        for (key in params.preset.do_merging.keySet()) {
            if ((key in params.do_merging.keySet()) && (params.preset.do_merging[key] instanceof Map)) {
                params.preset.do_merging[key].keySet().findAll {
                    !(it in params.do_merging[key].keySet())
                }.each {
                    presetsDoMerging[key] = presetsDoMerging.get(key, [:]); presetsDoMerging[key][it] = params.preset.do_merging[key][it]
                }
            } else if (!(key in params.do_merging.keySet())) {
                presetsDoMerging[key] = params.preset.do_merging[key]
            }
        }

        // add the presetsDoMerging to the params.do_merging, but if the preset is a map, merge it with the existing map
        for (key in presetsDoMerging.keySet()) {
            if (key in params.do_merging.keySet() && params.do_merging[key] instanceof Map) {
                params.do_merging[key].putAll(presetsDoMerging[key])
            } else {
                params.do_merging[key] = presetsDoMerging[key]
            }
        }
        return
    }

    // shared function for loading assay specific preset params
    private static <Type> void loadAssayPresetParams(params, assayParams, assayPresets, assay) {
        def presets = [:]
        for (key in assayPresets.keySet()) {
            if ((key in assayParams.keySet()) && (assayPresets[key] instanceof Map)) {
                assayPresets[key].keySet().findAll {
                    !(it in assayParams[key].keySet())
                }.each {
                    presets[key] = presets.get(key, [:]); presets[key][it] = assayPresets[key][it]
                }
            } else {
                if (!(key in assayParams.keySet())) {
                    presets[key] = assayPresets[key]
                }
            }
        }
        // add the presets to the assayParams, but if the preset is a map, merge it with the existing map
        for (key in presets.keySet()) {
            if (key in assayParams.keySet() && assayParams[key] instanceof Map) {
                assayParams[key].putAll(presets[key])
                params.core[key].putAll(presets[key])
                if ((assay == 'rna') || (assay == 'cytometry')) {
                    params.do_merging[key].putAll(presets[key])
                }
            } else {
                assayParams[key] = presets[key]
                params.core[key] = presets[key]
                if ((assay == 'rna') || (assay == 'cytometry')) {
                    params.do_merging[key] = presets[key]
                }
            }
        }
        return
    }

    private static <Type> void loadSampleLevelParams(assayParams) {
        // If no sample sheet was provided, we still need to associate sample-level params with the sample ID
        // Create a map to hold the sample overrides
        def samplesMap = [:]
        for (int i = 0; i < assayParams.sampleIds.size(); i++) {
            String sampleId = assayParams.sampleIds[i]
            def sampleLevelParams = ['barcode', 'cell', 'bead', 'gelBeadLambda',
                                     'beadMergeUmiThreshold', 'umiHamming', 'skipMerge',
                                     'recomFilter', 'trimFiveUmi', 'total_cDNA_reads',
                                     'trimThreePrime', 'trimFivePrime', 'mergeMethod', 'jaccard']
            for (key in sampleLevelParams) {
                // If the key is not already in the samplesMap, create it
                samplesMap[key] = samplesMap[key] ?: [:]
                if (assayParams[key] == null) {
                    // If the key is not set in the assayParams, skip it
                    if (key == "barcode") { // groovylint-disable-line
                        samplesMap[key]['force'] = samplesMap[key]['force'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['forceType'] = samplesMap[key]['forceType'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['filterType'] = samplesMap[key]['filterType'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['force'][sampleId] = null // groovylint-disable-line
                        samplesMap[key]['forceType'][sampleId] = null // groovylint-disable-line
                        samplesMap[key]['filterType'][sampleId] = null // groovylint-disable-line
                    } else if (key == "total_cDNA_reads") { // groovylint-disable-line
                        samplesMap[key][sampleId] = null
                    } else if (key == 'jaccard') { // groovylint-disable-line
                        samplesMap[key]['force'] = samplesMap[key]['force'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['force'][sampleId] = null // groovylint-disable-line
                    }
                    continue
                }
                // If the key is a map, iterate over its keys and add them to the samplesMap
                if (assayParams[key] instanceof LinkedHashMap) {
                    for (nestedKey in assayParams[key].keySet()) {
                        // If the nestedKey is not in the samplesMap[key], add it
                        samplesMap[key][nestedKey] = samplesMap[key][nestedKey] ?: [:]
                        // Add the sampleId to the nested map
                        samplesMap[key][nestedKey][sampleId] = assayParams[key][nestedKey]
                    }
                } else {
                    // If the key is not a map, just add it to the samplesMap
                    samplesMap[key][sampleId] = assayParams[key]
                }
            }
        }
        // put the samplesMap into assayParams
        assayParams.putAll(samplesMap)
    }

    private static <Type> void loadSampleSheetParams(assayParams) {
        // Update the assayParams object with the sample sheet metadata
        // Create a map to hold the sample overrides
        def samplesMap = [:]
        def keyToGrab
        for (int i = 0; i < assayParams.sampleSheetMetaMap.size(); i++) {
            def currentMap = assayParams.sampleSheetMetaMap[i][0]
            def ignore = [ 'sampleId', 'fastq1', 'fastq2'] // keys to ignore
            // iterate over keys of map, key/value is != [], store the key and value in assayParams
            // skip key if == sampleId, fastq1, or fastq2
            for (key in currentMap.keySet() - ignore) {
                // need to create a map for entry into assayParams
                String sampleId = currentMap['sampleId'] // groovylint-disable-line
                samplesMap[key] = samplesMap[key] ?: [:]
                if (currentMap[key] != []) {
                    // Sample override found
                    // Key to grab can be any type, set equal to currentMap[key]
                    keyToGrab = currentMap[key]
                } else {
                    // use the preset value if it exists
                    if ((assayParams[key]) == null && (key == "barcode")) { // groovylint-disable-line
                        samplesMap[key]['force'] = samplesMap[key]['force'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['force'][sampleId] = null // groovylint-disable-line
                        samplesMap[key]['forceType'] = samplesMap[key]['forceType'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['forceType'][sampleId] = null // groovylint-disable-line
                        samplesMap[key]['filterType'] = samplesMap[key]['filterType'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['filterType'][sampleId] = null // groovylint-disable-line
                        continue
                    } else if ((assayParams[key] == null) && (key == "cell")) { // groovylint-disable-line
                        samplesMap[key]['loaded'] = samplesMap[key]['loaded'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['loaded'][sampleId] = null // groovylint-disable-line
                        continue
                    } else if ((assayParams[key] == null) && (key == "bead")) { // groovylint-disable-line
                        samplesMap[key]['format'] = samplesMap[key]['format'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['length'] = samplesMap[key]['length'] ?: [:] // groovylint-disable-line
                        samplesMap[key]['format'][sampleId] = null // groovylint-disable-line
                        samplesMap[key]['length'][sampleId] = null // groovylint-disable-line
                        continue
                    }
                    else {
                        keyToGrab = assayParams[key]
                    }
                }
                // if not a nested parameter, add it to the samplesMap
                if (!(keyToGrab instanceof LinkedHashMap)) {
                    samplesMap[key][sampleId] = keyToGrab
                } else {
                    for (nestedKey in keyToGrab.keySet()) {
                        // if the nestedKey is not in the samplesMap[key], add it
                        samplesMap[key][nestedKey] = samplesMap[key][nestedKey] ?: [:]
                        if (keyToGrab[nestedKey] != []) {
                            // add the sampleId to the nested map
                            samplesMap[key][nestedKey][sampleId] = keyToGrab[nestedKey]
                        } else {
                            // Check if there is a preset value for the nested key
                            samplesMap[key][nestedKey][sampleId] = assayParams[key][nestedKey] ?: null
                        }
                    }
                }
            }
        }

        // put the samplesMap into assayParams
        assayParams.putAll(samplesMap)
    }

    // Function for extracting prefixes from input read files
    private static String getPrefix(path, pattern, log) {
        // Formatting input variables
        File file = new File(path)
        String name = file.name
        String glob = pattern

        // Finding location of glob wildcards
        int indexWildcard = glob.findIndexOf { it == '*' || it == '?' }
        int indexBracket = glob.findIndexOf { it == '{' || it == '[' }

        // Processing globs without wildcards
        if (indexWildcard == -1 && indexBracket == -1) {
            // Simplify file name if file matches glob exactly
            if (name == glob) {
                return file.simpleName
            // Error if glob doesn't match file name or have wildcards
            }
            log.error("ERROR: Not a valid file glob pattern: pattern = $glob file = $name")
            System.exit(1)
        }

        // Count wildcards before parameter expansion characters(brackets)
        int groups = glob.substring(0, indexBracket).findIndexValues { it == '*' || it == '?' } .size()

        // Convert to java glob
        def regex = glob // groovylint-disable-line
                .replace('.', '\\.')
                .replace('*', '(.*)')
                .replace('?', '(.?)')
                .replace('{', '(?:') // groovylint-disable-line
                .replace('}', ')')
                .replace(',','|') // groovylint-disable-line

        // Define pattern matcher
        def matcher = (name =~ /$regex/)

        // Initialize matcher and extract prefix from files that match glob
        if (matcher.matches()) {
            int end = matcher.end(groups)
            String prefix = name.substring(0, end)
            if (prefix.endsWith('-') || prefix.endsWith('_') || prefix.endsWith('.')) {
                prefix = prefix[0..-2]
            }

            return prefix
        }
        return null
    }

}
