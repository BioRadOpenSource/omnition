class Core {

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
            nextflow run digitalbiology/omnition -params-file <path to parameters file>

            Nextflow arguments:
            -params-file            Path to the parameter definition file
            -profile                Execute using Docker containers or Singularity containers(default: Singularity).
            -resume                 Resume a previous pipeline execution using the cached results(default: no resume).
            -r, -revision           Specify git tag, branch, or commit hash from which \
                to execute the pipeline(default: develop).

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

    // Function(generic) to format params accepting one or more values for input to channels
    static <Type> Type paramFormat(param) {
        if (param.class == String) {
            return param
        }
        return param.collect { it.value }
    }

    // Function for counting the number of elements in an object
    static int countElement(object) {
        int value

        if (object.class == String) {
            value = 1
        } else if (object.class == ArrayList) {
            value = object.size()
        }

        return value
    }

    // Function to check output s3 path for Tower runs
    static <Type> Type validateOutputS3(params, assayParams, profile, log) {
        String value
        Type path
        // Ignore validation for demo data that are smoke tests
        if (assayParams.output =~ /(smoke_tests)/) {
            log.info("[$assayParams.assay] Running smoke tests with Seqera Platform, output path set to "
                + assayParams.output + ".")
            return
        }

        if (profile =~ /(demo)/) {
            // Ignore validation for demo data on Seqera Platform
            String s3UploadPath
            if (params.options.output != null) {
                s3UploadPath = params.options.output
            }
            log.info("[$assayParams.assay] Using demo data with Seqera Platform, output path set to "
                + s3UploadPath + ".")
            assayParams.output = s3UploadPath
            return
        }

        // Condition to ignore qcp analysis path validation
        if (assayParams.output =~ /(qcp_omnition_analysis_results)/) {
            // Ignore validation for qcp data on Seqera Platform
            String s3UploadPath
            if (params.options.output != null) {
                s3UploadPath = params.options.output
            }
            log.info("[$assayParams.assay] Analzying QCP data with Seqera Platform, output path set to "
                + s3UploadPath + ".")
            assayParams.output = s3UploadPath
            return
        }

        // Check if output set correctly
        if (assayParams.output != null) {
            value = assayParams.output
            path = value.split(params.tower.aws_analysis_bucket)
            // condition 1: if output path is exactly the minimum path
            // condition 2: if output path contains anything but the minimum path at beginning
            if ((path.size() == 0) || (path[0] != "")) { // groovylint-disable-line
                log.error("[$assayParams.assay] Invalid output path: $value")
                System.exit(1)
            }
        }
    }

    // Function to check 'output' parameter
    static String validateOutput(params, assayParams, log) {
        String value
        File file

        // Check if a specific output directory is set
        if (assayParams.output != null) {
            value = assayParams.output
            if (value.endsWith("/")) {
                value = value.substring(0, value.length() - 1) // remove trailing '/'
            }
            file = new File(value)
        } else {
            value = params.options.output // this is the default "./results" dir
            file = new File(value)
        }
        // Check if the specified directory exists and if has anything in it other than the trace file
        if ((file.list() != null) && (file.list().length > 0) && !onlyTraceFileInOutDir(file, params)) {
            // Check if the force flag was not supplied
            if (params.force == null) {
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
    static <Type> Type validateParamsFile(params, assayMode) {
        Type value

        // If there are assay level params specified in the CLI, grab them
        if (params.options.paramsFile?.get(assayMode.toLowerCase())) {
            value = params.options.paramsFile
        } else {
            // If assay-level parameters are not found from the CLI params
            // then it means no params-file was given
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

    // Function to check the 'mixed' parameter
    static boolean validateMixed(params, log, messages) {
        boolean value

        // Check if set
        if (params.mixed != null) {
            // Check if class type is correct
            if (params.mixed.class != Boolean) {
                log.error("ERROR: [$params.assay] The 'mixed' parameter must be boolean(true or false).")
                System.exit(1)
            } else {
                log.info("INFO: [$params.assay] Executing mixed species workflow.")
                messages.add("INFO: [$params.assay] Executing mixed species workflow.")
                value = params.mixed
            }
        } else {
            value = false
        }
        return value
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

    // Function to check the reference 'fasta' file parameter
    static <Type> Type validateReferenceFasta(params, log) {
        Type value
        int count

        // Check if set
        if (params.reference?.fasta != null) {
            value = Core.paramFormat(params.reference?.fasta)
            // Check if endings are valid
            if (value.stream().allMatch { it.matches('\\S+(.fa|.fa.gz|.fna|.fna.gz|.fasta|.fasta.gz)') }) { // groovylint-disable-line
                // Check if files exist
                if (value.stream().allMatch { new File(it).exists() }) {
                    count = Core.countElement(value)
                    // Check if valid number of files have been provided
                    if (params.mixed != true && count != 1) {
                        log.error("ERROR: [$params.assay] The reference 'fasta' parameter requires one file.\
                            If using \
                            two references, set 'mixed: true' in the parameters file.")
                        System.exit(1)
                    } else if (params.mixed == true && count != 2) {
                        log.error("ERROR: [$params.assay] The reference 'fasta' parameter requires two files when \
                            'mixed: true' is set in the paramters file.")
                        System.exit(1)
                    } else {
                        return value
                    }
                } else {
                    log.error("ERROR: [$params.assay] The reference'fasta' parameter file does not exist: $value")
                    System.exit(1)
                }
            } else {
                log.error("ERROR: [$params.assay] The reference 'fasta' parameters must end in '.fa', \
                    '.fa.gz', '.fna', '.fna.gz', '.fasta', or '.fasta.gz'.")
                System.exit(1)
            }
        } else {
            log.error("ERROR: [$params.assay] Must set the reference 'fasta' parameter.")
            System.exit(1)
        }
    }

    // Function to check the reference 'gtf' file parameter
    static <Type> Type validateReferenceGtf(params, log) {
        Type value
        int count

        // Check if set
        if (params.reference?.gtf != null) {
            value = Core.paramFormat(params.reference?.gtf)
            // Check if endings are valid
            if (value.stream().allMatch { it.matches('\\S+(.gtf|.gtf.gz)') }) {
                // Check if files exist
                if (value.stream().allMatch { new File(it).exists() }) {
                    count = Core.countElement(value)
                    // Check if valid number of files have been provided
                    if (params.mixed != true && count != 1) {
                        log.error("ERROR: [$params.assay] The reference 'gtf' parameter requires one file. \
                            If using two references, set 'mixed: true' in the parameters file.")
                        System.exit(1)
                    } else if (params.mixed == true && count != 2) {
                        log.error("ERROR: [$params.assay] The reference 'gtf' parameter requires two \
                            files when 'mixed: true' is set in the paramters file.")
                        System.exit(1)
                    } else {
                        return value
                    }
                } else {
                    log.error("ERROR: [$params.assay] The reference 'gtf' parameter file does not exist: $value")
                    System.exit(1)
        }
            } else {
                log.error("ERROR: [$params.assay] The reference 'gtf' parameters must end in '.gtf' or '.gtf.gz'.")
                System.exit(1)
    }
        } else {
            log.error("ERROR: [$params.assay] Must set the reference 'gtf' parameter.")
            System.exit(1)
        }
    }

    // Function for setting the overrides for specific samples if needed
    static boolean validateSampleOverride(params, log) {
        // Initializing variables
        boolean value
        // Iterating over the provided override sample IDs
        if (params.overrides != null) {
            Set overridesamples = params.overrides.keySet()

            for (int i = 0; i < overridesamples.size(); i++) {
                String key = overridesamples[i]
                // Check if class type is correct
                if (key.class != String) {
                    log.error("ERROR: [$params.assay] The 'overrides' parameter \
                        $key must be a String.")
                    System.exit(1)
                } else if (params.sampleIds.contains(key)) { // Valid Sample
                    value = true
                } else {
                    log.error("ERROR: [$params.assay] The 'overrides' parameter \
                        $key is not in the sample list: $params.sampleIds.")
                    System.exit(1)
                }
            }
        } else {
            value = false
        }
        return value
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
        Core.paramError(params.assay, 'workflow', log)
        return
    }

    // Function to ensure the input file names match the Illumina naming convention
    // groovylint-disable
    //https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm 
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
                    it.toString().matches('.*\\/Undetermined[\\w\\.\\-]+') } //groovylint-disable-line
            // Warn if FASTQ files beginning with 'Undetermined' are found
            if (test) {
                log.warn("[$params.assay] Undetermined FASTQ read files will be ignored.")
                messages.add("WARN: [$params.assay] Undetermined FASTQ read files will be ignored.")
            }
            if (sampleFiles.size() != 0) {
                if (!Core.validateInputNamesAWS(sampleFiles)) {
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
            Core.paramError(params.assay, 'input', log)
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
                if (!Core.validateInputNames(file)) {
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
            Core.paramError(params.assay, 'input', log)             
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

    static List<String> getSpecies(assay, gtf, log) {
        // Creating list of species names for all files matching all globs
        List speciesList = gtf
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
                    + "parameters and file name requirements.")
            System.exit(1)
        } else {
            return speciesList
        }
    }

    // Function for extracting all sample IDs matching a given pattern(NOTE: filters out "Undetermined" samples)
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

        // Creating list of sample IDs for all files matching  all globs
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

    // AWS Batch implementation
    static List<String> getSampleIdsAWS(assay, sampleFiles, pattern, log, messages) {
        List globList = pattern instanceof List ? pattern : [pattern]

        // Creating list of all file prefixes
        List prefixList = []
        for (int i = 0; i < globList.size(); i++) {
            List matchList = sampleFiles
                .stream()
                .map(j -> getPrefix(j.toString(), globList[i].toString(), log))
                .collect()

            prefixList.addAll(matchList)
        }

        // Creating list of sample IDs for all files matching  all globs
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

    // Function for checking if the only file in the output directory is a trace file
    private static boolean onlyTraceFileInOutDir(dir, params) {
        // Gather information on the output dir files and count differences between
        // expected files
        List expected = [".nextflow", ".nextflow.log", "pipeline_info", "work"]
        int diff = (dir.list() - expected).size()

        // Get the pipeline info subDir to check for only trace file existing
        int index = dir.list().findIndexValues { it ==~ "pipeline_info" } [0] // groovylint-disable-line
        File subDir = dir.listFiles()[index]

        // Get trace name string that allows for setting prefix
        String trace = params.prefix ? (params.prefix + "-omnition-trace-") : ("omnition-trace-")

        if (diff > 0 || !subDir.directory) {
            return false
        }
        if (subDir.list().length > 1) {
            return false
        }
        return subDir.listFiles()[0].toString().contains(trace)
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
