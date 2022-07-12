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

    // Function to create help message
    static String helpMessage(monochrome_logs) { // groovylint-disable-line
        // logo = bioradLogo(monochrome_logs)

        def message = '' // groovylint-disable-line

        message = String.format(
            '''\n
            Standard usage
            nextflow run BioRadOpenSource/omnition-test -params-file <path to parameters file>

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

    // Function to check 'outputDir' parameter
    static String validateOutputDir(params, log) {
        String value
        File file

        // Check if a specific outputDir directory is set
        if (params.outputDir) {
            value = params.outputDir
            file = new File(value)
            // Check if the specified directory exists and if has anything in it
            if ((file.list() != null) && (file.list().length > 0)) {
                // Check if the force flag was not supplied
                if (!params.force) {
                    log.error("ERROR: The specified output directory is not empty. \
                        This can be overridden with the --force flag: $params.outputDir.")
                    System.exit(1)
                }
            }
        }
    }

    // Function to ensure the input file names match the Illumina naming convention
    // groovylint-disable
    //https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/NamingConvention_FASTQ-files-swBS.htm 
    // groovylint-enable
    static boolean validateInputNames(directory) {
        return directory.listFiles().stream().allMatch {
                    it.name.matches('[a-zA-Z0-9]+_S[0-9]+_L[0-9]{3}_R[0-9]+_001.f(ast)?q.gz') ||
                     (!it.name.endsWith('fastq.gz') &&
                     !it.name.endsWith('fq.gz')) ||
                     it.directory }
    }

    // Function to ensure that species name / prefix match between fasta-gtf pairs
    static void validateReferenceNames(fasta, gtf, log) {
        // Get just the prefix from the file names
        List gtfSpeciesList = gtf
            .stream()
            .filter(Objects::nonNull)
            .map(i -> i - ~/^\/.*\//)
            .map(i -> i - ~/\..*$/)
            .distinct()
            .collect()

        List fastaSpeciesList = fasta
            .stream()
            .filter(Objects::nonNull)
            .map(i -> i - ~/^\/.*\//)
            .map(i -> i - ~/\..*$/)
            .distinct()
            .collect()

        // The prefixes should match and be in the same order
        if (gtfSpeciesList == (fastaSpeciesList)) {
            log.info('INFO: fasta and gtf match')
        } else {
            log.error('ERROR: fasta and gtf prefixes do not match, exiting')
            System.exit(1)
        }
    }

    static List<String> getSpecies(assay, gtf, log) {
        // Creating list of species names for all files matching all globs
        List speciesList = gtf
            .stream()
            .filter(Objects::nonNull)
            .map(i -> i - ~/^\/.*\//)
            .map(i -> i - ~/\..*$/)
            .distinct()
            .collect()

        // Checking if files were found
        if (speciesList.empty) {
            log.error("ERROR: [$assay] No GTF files found in input directory. Check \
                parameters and file name requirements.")
            System.exit(1)
        } else {
            return speciesList
        }
    }

    // Function for extracting all sample IDs matching a given pattern(NOTE: filters out "Undetermined" samples)
    static List<String> getSampleIds(assay, directory, pattern, log) {
        // Formatting inputs
        List fileList = new File(directory).listFiles()
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
            log.error("ERROR: [$assay] No FASTQ files found in input directory. \
                Check parameters and file name requirements.")
            System.exit(1)
        } else {
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
                .replace(',','|')

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
