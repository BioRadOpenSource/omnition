/* groovylint-disable */
import Core
/* groovylint-enable */

class Atac {

    //Function to validate AWS reference directory for AWS
    static String validateReferenceDirectoryAWS(params, log) {
        String value

        // Check if set
        if (params.atac.reference?.directory) {
            value = params.atac.reference.directory
            return value
        }
        log.error("ERROR: [$params.atac.assay] Must set the reference 'directory' parameter.")
        System.exit(1)
        return
    }

    // Function to check the contaminant 'directory' parameter for AWS
    static <Type> Type validateContaminantDirectoryAWS(params, log) {
        Type value
        log

        // Check if set
        if (params.atac.contaminant?.directory) {
            value = params.atac.contaminant.directory
            return value
        }
        return false
    }

    // Function to check the 'workflow' parameter
    static String validateWorkflow(params, log) {
        String value

        // Check if set
        if (params.atac.workflow) {
            value = params.atac.workflow
            // Check if provided value is valid
            if (!(value in [ 'reference', 'analysis', 'full' ])) {
                log.error("ERROR: [$params.atac.assay] The 'workflow' parameter must be \
                    'reference', 'analysis', or 'full'.")
                System.exit(1)
            }
            return value
        }
        Core.paramError(params.atac.assay, 'workflow', log)
        return
    }

    // Function to check the reference 'directory' parameter
    static String validateReferenceDirectory(params, log) {
        String value

        // Check if set
        if (params.atac.reference?.directory) {
            value = params.atac.reference.directory
            return value
        }
        log.error("ERROR: [$params.atac.assay] Must set the reference 'directory' parameter.")
        System.exit(1)
        return
    }

    // Function to check the reference 'fasta' file parameter
    static <Type> Type validateReferenceFasta(params, log) {
        Type value
        int count

        // Check if set
        if (params.atac.reference?.fasta) {
            value = Core.paramFormat(params.atac.reference?.fasta)
            // Check if endings are valid
            if (value.stream().allMatch { it.matches('\\S+(.fa|.fa.gz|.fna|.fna.gz|.fasta|.fasta.gz)') }) {
                // Check if files exist
                if (value.stream().allMatch { new File(it).exists() }) {
                    count = Core.countElement(value)
                    // Check if valid number of files have been provided
                    if (!params.atac.mixed && count != 1) {
                        log.error("ERROR: [$params.atac.assay] The reference 'fasta' parameter requires one file.\
                            If using \
                            two references, set 'mixed: true' in the parameters file.")
                        System.exit(1)
                    } else if (params.atac.mixed && count != 2) {
                        log.error("ERROR: [$params.atac.assay] The reference 'fasta' parameter requires two files when \
                            'mixed: true' is set in the paramters file.")
                        System.exit(1)
                    } else {
                        return value
                    }
                } else {
                    log.error("ERROR: [$params.atac.assay] The reference'fasta' parameter file does not exist: $value")
                    System.exit(1)
        }
            } else {
                log.error("ERROR: [$params.atac.assay] The reference 'fasta' parameters must end in '.fa', \
                    '.fa.gz', '.fna', '.fna.gz', '.fasta', or '.fasta.gz'.")
                System.exit(1)
    }
        } else {
            log.error("ERROR: [$params.atac.assay] Must set the reference 'fasta' parameter.")
            System.exit(1)
}
    }

    // Function to check the reference 'gtf' file parameter
    static <Type> Type validateReferenceGtf(params, log) {
        Type value
        int count

        // Check if set
        if (params.atac.reference?.gtf) {
            value = Core.paramFormat(params.atac.reference?.gtf)
            // Check if endings are valid
            if (value.stream().allMatch { it.matches('\\S+(.gtf|.gtf.gz)') }) {
                // Check if files exist
                if (value.stream().allMatch { new File(it).exists() }) {
                    count = Core.countElement(value)
                    // Check if valid number of files have been provided
                    if (!params.atac.mixed && count != 1) {
                        log.error("ERROR: [$params.atac.assay] The reference 'gtf' parameter requires one file. \
                            If using two references, set 'mixed: true' in the parameters file.")
                        System.exit(1)
                    } else if (params.atac.mixed && count != 2) {
                        log.error("ERROR: [$params.atac.assay] The reference 'gtf' parameter requires two \
                            files when 'mixed: true' is set in the paramters file.")
                        System.exit(1)
                    } else {
                        return value
                    }
                } else {
                    log.error("ERROR: [$params.atac.assay] The reference 'gtf' parameter file does not exist: $value")
                    System.exit(1)
        }
            } else {
                log.error("ERROR: [$params.atac.assay] The reference 'gtf' parameters must end in '.gtf' or '.gtf.gz'.")
                System.exit(1)
    }
        } else {
            log.error("ERROR: [$params.atac.assay] Must set the reference 'gtf' parameter.")
            System.exit(1)
        }
    }

    // Function to check the contaminant 'directory' parameter
    static <Type> Type validateContaminantDirectory(params, log) {
        Type value
        log

        // Check if set
        if (params.atac.contaminant?.directory) {
            value = params.atac.contaminant.directory
            return value
        }
        return false
    }

    // Function to check the contaminant 'fasta' file parameter
    static <Type> Type validateContaminantFasta(params, log) {
        Type value

        // Check if set
        if (params.atac.contaminant?.fasta) {
            value = Core.paramFormat(params.atac.contaminant?.fasta)
            // Check if endings are valid
            if (value.stream().allMatch { it.matches('\\S+(.fa|.fa.gz|.fna|.fna.gz|.fasta|.fasta.gz)') }) { // groovylint-disable-line
                // Check if files exist
                if (value.stream().allMatch { new File(it).exists() }) {
                    // Check if contaminant 'directory' has been set
                    if (params.atac.contaminant?.directory) {
                        return value
                    }
                    log.error("ERROR: [$params.atac.assay] The contaminant 'directory' parameter must \
                        be set if providing contaminants.")
                    System.exit(1)
                } else {
                    log.error("ERROR: [$params.atac.assay] The contaminant 'fasta' parameter file \
                        does not exist: $value")
                    System.exit(1)
        }
            } else {
                log.error("ERROR: [$params.atac.assay] The contaminant 'fasta' parameters must end in \
                    '.fa', '.fa.gz', '.fna', '.fna.gz', '.fasta', or '.fasta.gz'.")
                System.exit(1)
    }
        } else {
            return false
        }
    }

    // Function to check the 'mixed' parameter
    static boolean validateMixed(params, log, messages) {
        boolean value

        // Check if set
        if (params.atac.mixed) {
            // Check if class type is correct
            if (params.atac.mixed.class != Boolean) {
                log.error("ERROR: [$params.atac.assay] The 'mixed' parameter must be boolean(true or false).")
                System.exit(1)
            } else {
                log.info("INFO: [$params.atac.assay] Executing mixed species workflow.")
                messages.add("INFO: [$params.atac.assay] Executing mixed species workflow.")
                value = params.atac.mixed
            }
        } else {
            value = false
        }
        return value
    }

    // Function to check the 'tierroroverride' parameter
    static boolean validateTIErrorOverride(params, log, messages) {
        boolean value

        // Check if set
        if (params.atac.tierroroverride) {
            // Check if class type is correct
            if (params.atac.tierroroverride.class != Boolean) {
                log.error("ERROR: [$params.atac.assay] The 'tierroroverride' parameter must be boolean(true or false).")
                System.exit(1)
            } else {
                log.warn("[$params.atac.assay] Overriding errors in the config file. Only TIs specified \
                    in the config file will be used.")
                messages.add("WARN: [$params.atac.assay] Overriding errors in the config file. \
                    Only TIs specified in the config file will be used.")
                value = params.atac.tierroroverride
            }
        } else {
            value = params.preset.atac.tierroroverride
        }
        return value
    }

    // Function to check the 'barcodedTn5' parameter
    static boolean validateBarcodedTn5(params, log, messages) {
        boolean value

        // Check if set
        if (params.atac.barcodedTn5) {
            // Check if class type is correct
            if (params.atac.barcodedTn5.class != Boolean) {
                log.error("ERROR: [$params.atac.assay] The 'barcodedTn5' parameter must be boolean(true or false).")
                System.exit(1)
            } else {
                log.info("INFO: [$params.atac.assay] Executing barcodedTn5 workflow.")
                messages.add("INFO: [$params.atac.assay] Executing barcodedTn5 workflow.")
                value = params.atac.barcodedTn5
            }
        } else {
            value = false
        }
        return value
    }

    // Function to check the barcodedTn5 config parameter
    static <Type> Type validateBarcodedTn5Config(params, log, messages) {
        Type value

        // Check if barcodedTn5 set
        if (params.atac.barcodedTn5) {
            // Check if parameter set
            if (params.atac.barcodedTn5Config) {
                value = params.atac.barcodedTn5Config
                // Check if file exists
                if (new File(value).exists()) {
                    // Check if file is a csv
                    if (value.stream().allMatch { it.matches('\\S+(.csv)') }) {
                        return value
                    }
                    log.error("ERROR: [$params.atac.assay] The 'BarcodedTn5Config' parameter file \
                        must be a csv file $value")
                    System.exit(1)
                } else {
                    log.error("ERROR: [$params.atac.assay] The 'BarcodedTn5Config' parameter \
                        file does not exist: $value")
                    System.exit(1)
            }
            } else {
                log.info("INFO: [$params.atac.assay] No 'BarcodedTn5Config' Config file provided, \
                    executing superloading workflow.")
                messages.add("INFO: [$params.atac.assay] No 'BarcodedTn5Config' Config file provided, \
                    executing superloading workflow.")
                return false
        }
        } else {
            return false
    }
    }

    // Function to check the TI parameter
    static <Type> Type validateTI(params, log, messages) {
        Type value

        // Check if barcodedTn5 set
        if (params.atac.barcodedTn5) {
            // Check if parameter set
            if (params.atac.ti) {
                value = params.atac.ti
            } else {
                log.info("INFO: [$params.atac.assay] No TI list provided, using preset TIs.")
                messages.add("INFO: [$params.atac.assay] No TI list provided, using preset TIs.")
                value = params.preset.atac.ti
            }

            return(value)
        }
        return false
    }

    // Function to check the TIread parameter
    static <Type> Type validateTIRead(params, log, messages) {
        Type value

        // Check if barcodedTn5 set
        if (params.atac.barcodedTn5) {
            // Check if parameter set
            if (params.atac.tiread) {
                value = params.atac.tiread

                // Check if i7asti is false
                if (params.atac.i7asti) {
                    log.error("ERROR: [$params.atac.assay] The 'tiread' parameter cannot be used in \
                        conjuction with the 'i7asti' parameter.")
                }

                // Performing type checks
                if (!value.matches('(r1|r2)')) {
                    log.error("ERROR: [$params.atac.assay] The 'tiread' parameter must be 'r1', or 'r2'.")
                }
            } else {
                value = params.preset.atac.tiread
            }

            if (params.atac.i7asti) {
                log.info("INFO: [$params.atac.assay] TIs are set to be located in the i7 barcode.")
                messages.add("INFO: [$params.atac.assay] TIs are set to be located in the i7 barcode.")
            } else {
                log.info("INFO: [$params.atac.assay] TIs are set to be located in $value.")
                messages.add("INFO: [$params.atac.assay] TIs are set to be located in $value.")
            }
            return(value)
        }
        return false
    }

    // Function to check the 'i7asti' parameter
    static boolean validatei7asti(params, log, messages) {
        boolean value

        // Check if set
        if (params.atac.i7asti) {
            // Check if class type is correct
            if (params.atac.i7asti.class != Boolean) {
                log.error("ERROR: [$params.atac.assay] The 'i7asti' parameter must be boolean(true or false).")
                System.exit(1)
            } else if (params.atac.barcodedTn5 != true) {
                log.error("ERROR: [$params.atac.assay] The 'BarcodedTn5' parameter must be true \
                    in order to set 'i7asti' parameter to true.")
                System.exit(1)
            } else {
                String message = "INFO: [$params.atac.assay] Using i7 sequence as a tagmentation index."
                log.info(message)
                messages.add(message)
                value = params.atac.i7asti
            }
        } else {
            value = false
        }
        return value
    }

    // Function to check the 'config' parameter
    static String validateConfig(params, log) {
        String value

        // Check if set
        if (params.atac.config) {
            value = params.atac.config
            // Check if file exists
            if (new File(value).exists()) {
                return value
            }
            log.error("ERROR: [$params.atac.assay] The 'config' parameter file does not exist: $value")
            System.exit(1)
        } else {
            Core.paramError(params.atac.assay, 'config', log)
        }
    }

    // Function for setting the overrides for specific samples if needed
    static boolean validateSampleOverride(params, log) {
        // Initializing variables
        boolean value
        // Iterating over the provided override sample IDs
        if (params.atac.overrides) {
            Set overridesamples = params.atac.overrides.keySet()

            for (int i = 0; i < overridesamples.size(); i++) {
                String key = overridesamples[i]
                // Check if class type is correct
                if (key.class != String) {
                    log.error("ERROR: [$params.atac.assay] The 'overrides' parameter \
                        $key must be a String.")
                    System.exit(1)
                } else if (params.atac.sampleIds.contains(key)) { // Valid Sample
                    value = true
                } else {
                    log.error("ERROR: [$params.atac.assay] The 'overrides' parameter \
                        $key is not in the sample list: $params.atac.sampleIds.")
                    System.exit(1)
                }
            }
        } else {
            value = false
        }
        return value
    }

    // Function for setting the overrides for specific TIs if needed
    static boolean validateTIOverride(params, log) {
        // Initializing variables
        boolean value
        // Iterating over the provided override sample IDs
        if (params.atac.overrides && params.atac.barcodedTn5) {
            Set overridesamples = params.atac.overrides.keySet()

            for (int i = 0; i < overridesamples.size(); i++) {
                String key = overridesamples[i]
                Set tinames = params.atac.ti.keySet()
                // Get parameter overrides for the sample
                def entries = params.atac.overrides.get(key)
                .collectEntries { entry -> [ entry.key, entry.value ] }
                for (e in entries) {
                    // If the override is a TI, class will be a Map not a string, where the key is the TI
                    if (e.value.getClass() != String) {
                        if (tinames.contains(e.key)) { // Valid TI
                            value = true
                        } else {
                            log.error("ERROR: [$params.atac.assay] The 'overrides' parameter \
                                $e.key is not in the ti list: $tinames.")
                            System.exit(1)
                        }
                    }
                }
            }
        } else {
            value = false
        }
        return value
    }

    // Function to check 'input' directory parameter
    static String validateInput(params, log, messages) {
        String value
        File file

        // Check if set
        if (params.atac.input) {
            value = params.atac.input
            file = new File(value)
            // Check if directory exists
            if (file.exists()) {
                boolean test = file.listFiles().stream().anyMatch {
                    it.toString().matches('.*\\/Undetermined[\\w\\.\\-]+') }
                // Warn if FASTQ files beginning with 'Undetermined' are found
                if (test) {
                    log.warn("[$params.atac.assay] Undetermined FASTQ read files will be ignored.")
                    messages.add("WARN: [$params.atac.assay] Undetermined FASTQ read files will be ignored.")
                }
                if (!Core.validateInputNames(file)) {
                    log.error("ERROR: [$params.atac.assay] All FASTQ files in the 'input' parameter directory must" +
                    " match the Illumina naming convention SampleName_S1_L001_R1_001.fastq.gz")
                    System.exit(1)
                }
                return value
            }
            log.error("ERROR: [$params.atac.assay] The 'input' parameter directory does not exist: $value")
            System.exit(1)
        } else {
            Core.paramError(params.atac.assay, 'input', log)
        }
    }

    // Function to check 'contaminant' file parameter
    static <Type> Type validateReferenceBlocklist(params, log) {
        Type value
        int count

        // Check if set
        if (params.atac.reference?.blocklist) {
            value = Core.paramFormat(params.atac.reference?.blocklist)
            // Check if endings are valid
            if (value.stream().allMatch { it.matches('\\S+(blocklist.bed)') }) {
                if (value.stream().allMatch { new File(it).exists() }) {
                    count = Core.countElement(value)
                    // Check if valid number of files have been provided
                    if (!params.atac.mixed && count != 1) {
                        log.error("ERROR: [$params.atac.assay] The reference 'blocklist' parameter accepts at most \
                            one file per reference genome. If using two references, set \
                            'mixed: true' in the parameters file.")
                        System.exit(1)
                    } else if (params.atac.mixed && count > 2) {
                        log.error("ERROR: [$params.atac.assay] The reference 'blocklist' parameter accepts \
                            at most one file per reference genome.")
                        System.exit(1)
                    } else {
                        return value
                    }
                } else {
                    log.error("ERROR: [$params.atac.assay] The reference 'blocklist' parameter \
                        file does not exist: $value")
                    System.exit(1)
        }
            } else {
                log.error("ERROR: [$params.atac.assay] The reference 'blocklist' \
                    parameters must end in 'blocklist.bed'.")
                System.exit(1)
    }
        } else {
            // Check if workflow is analysis
            if (params.atac.workflow == 'analysis') {
                Core.paramError(params.atac.assay + '_' + params.atac.workflow, 'blocklist', log)
            } else {
                return false
            }
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
                .map(j -> Core.getPrefix(j.toString(), globList[i].toString(), log))
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
            log.error("ERROR: [$assay] No FASTQ files found in input directory. Check \
                parameters and file name requirements.")
            System.exit(1)
        } else {
            return sampleIdList
        }
    }

    // Function for setting cell settings and overriding if needed
    // ATAC does not currently use cells loaded, but keeping it in incase its used in the future
    static <Type> Type validateCell(params, log) {
        // Initializing variables
        List cellSettings = [ 'loaded' ]
        def cellMap = [:]

        // Iterating over setting types
        for (int i = 0; i < cellSettings.size(); i++) {
            String setting = cellSettings[i]
            def sampleMap = [:]

            // Iterating over sample IDs
            for (int j = 0; j < params.atac.sampleIds.size(); j++) {
                String key = params.atac.sampleIds[j]
                Type value

                if (params.atac.overrides?.get(key)?.cell?.get(setting)) {
                    value = params.atac.overrides.get(key).cell.get(setting)
                } else if (params.atac.cell?.get(setting)) {
                    value = params.atac.cell.get(setting)
                } else {
                    value = params.cell.get(setting)
                }

                // Performing type checks
                if (value != null && value.class != Integer) {
                    log.error("ERROR: [$params.atac.assay] Cell $setting parameter must \
                        be provided as an integer: $value")
                    System.exit(1)
                }

                // Add <key, value> pairs to map
                sampleMap.put(key, value)
            }

            // Add <key, value> pairs to map
            cellMap.put(setting, sampleMap)
        }

        return cellMap
    }

    // Function for setting barcode settings and overriding if needed
    static <Type> Type validateBarcode(params, log) {
        // Initializing variables
        List barcodeSettings = [ 'force' ]
        def barcodeMap = [:]

        // Is the run barcodedTn5?
        if (params.atac.barcodedTn5) {
            Set tinames = params.atac.ti.keySet()

            // Iterating over setting types
            for (int i = 0; i < barcodeSettings.size(); i++) {
                String setting = barcodeSettings[i]
                def sampleMap = [:]

                // Iterating over sample IDs
                for (int j = 0; j < params.atac.sampleIds.size(); j++) {
                    String key = params.atac.sampleIds[j]

                    // Iterating over tagmentation indexes
                    for (int k = 0; k < tinames.size(); k++) {
                        String tikey = tinames[k]
                        Type value

                        if (params.atac.overrides?.get(key)?.get(tikey)?.barcode?.get(setting)) {
                            value = params.atac.overrides.get(key).get(tikey).barcode.get(setting)
                        } else if (params.atac.overrides?.get(key)?.barcode?.get(setting)) {
                            value = params.atac.overrides.get(key).barcode.get(setting)
                        } else if (params.atac.barcode?.get(setting)) {
                            value = params.atac.barcode.get(setting)
                        } else {
                            value = params.barcode.get(setting)
                        }

                        // Performing type checks
                        if (value != null && value.class != Integer) {
                            log.error("ERROR: [$params.atac.assay] Barcode $setting parameter must \
                                be provided as an integer: $value")
                            System.exit(1)
                        }

                        // Add <key, value> pairs to map
                        sampleMap.put(key + '-' + params.atac.ti.get(tikey), value)
                    }
                }

                // Add <key, value> pairs to map
                barcodeMap.put(setting, sampleMap)
            }
        } else {
            // Iterating over setting types
            for (int i = 0; i < barcodeSettings.size(); i++) {
                String setting = barcodeSettings[i]
                def sampleMap = [:]

                // Iterating over sample IDs
                for (int j = 0; j < params.atac.sampleIds.size(); j++) {
                    String key = params.atac.sampleIds[j]
                    Type value

                    if (params.atac.overrides?.get(key)?.barcode?.get(setting)) {
                        value = params.atac.overrides.get(key).barcode.get(setting)
                    } else if (params.atac.barcode?.get(setting)) {
                        value = params.atac.barcode.get(setting)
                    } else {
                        value = params.barcode.get(setting)
                    }

                    // Performing type checks
                    if (value != null && value.class != Integer) {
                        log.error("ERROR: [$params.atac.assay] Barcode $setting parameter must be \
                            provided as an integer: $value")
                        System.exit(1)
                    }

                    // Add <key, value> pairs to map
                    sampleMap.put(key, value)
                }

                // Add <key, value> pairs to map
                barcodeMap.put(setting, sampleMap)
            }
        }

        return barcodeMap
    }

    // Function for setting the number of bases to trim from R2 reads
    static <Type> Type validateTrim(params, log) {
        // Initializing variables
        def trimMap = [:]

        // Is the run barcodedTn5?
        if (params.atac.barcodedTn5) {
            Set tinames = params.atac.ti.keySet()

            // Iterating over sample IDs
            for (int j = 0; j < params.atac.sampleIds.size(); j++) {
                String key = params.atac.sampleIds[j]

                // Iterating over tagmentation indexes
                for (int k = 0; k < tinames.size(); k++) {
                    String tikey = tinames[k]
                    Type value
                    if (params.atac.overrides?.get(key)?.get(tikey)?.trim) {
                        value = params.atac.overrides.get(key).get(tikey).trim
                    } else if (params.atac.overrides?.get(key)?.trim) {
                        value = params.atac.overrides.get(key).trim
                    } else if (params.atac.trim) {
                        value = params.atac.trim
                    } else {
                        value = params.preset.atac.trim
                    }

                    // Performing type checks
                    if (value != null && value.class != Integer) {
                        log.error("ERROR: [$params.atac.assay] The 'trim' parameter must be \
                            provided as an integer: $value")
                        System.exit(1)
                    }

                    // Add <key, value> pairs to map
                    trimMap.put(key + '-' + params.atac.ti.get(tikey), value)
                }
            }
        } else {
            // Iterating over sample IDs
            for (int j = 0; j < params.atac.sampleIds.size(); j++) {
                String key = params.atac.sampleIds[j]
                Type value

                if (params.atac.overrides?.get(key)?.trim) {
                    value = params.atac.overrides.get(key).trim
                } else if (params.atac.trim) {
                    value = params.atac.trim
                } else {
                    value = params.preset.atac.trim
                }

                // Performing type checks
                if (value != null && value.class != Integer) {
                    log.error("ERROR: [$params.atac.assay] The 'trim' parameter must be provided as an integer: $value")
                    System.exit(1)
                }

                // Add <key, value> pairs to map
                trimMap.put(key, value)
            }
        }

        return trimMap
    }

    // Function to check the 'mitoContig' parameter
    static <Type> Type validateMitoContig(params, log) {
        String value

        // Check if set
        if (params.atac.mitoContig) {
            // Check if class type is correct
            if (params.atac.mitoContig.class != String) {
                log.error("ERROR: [$params.atac.assay] The 'mitoContig' parameter must be a string.")
                System.exit(1)
            } else {
                value = params.atac.mitoContig
            }
        } else {
            value = params.preset.atac.mitoContig
        }
        return value
    }

    // Function to check the 'tssWindowSize' parameter
    static <Type> Type validateTSSWindow(params, log) {
        int value

        // Check if set
        if (params.atac.tssWindowSize) {
            // Check if class type is correct
            if (params.atac.tssWindowSize.class != Integer) {
                log.error("ERROR: [$params.atac.assay] The 'tssWindowSize' parameter must be an even integer.")
                System.exit(1)
            } else {
                // Check if given value is even
                if (params.atac.tssWindowSize % 2 == 0) {
                    value = params.atac.tssWindowSize
                } else {
                    log.error("ERROR: [$params.atac.assay] The 'tssWindowSize' parameter must be an even integer.")
                }
            }
        } else {
            value = params.preset.atac.tssWindowSize
        }
        return value
    }

    // Function to check the 'qualityThreshold' parameter
    static <Type> Type validateQualityThreshold(params, log) {
        int value

        // Check if set
        if (params.atac.qualityThreshold) {
            // Check if class type is correct
            if (params.atac.qualityThreshold.class != Integer) {
                log.error("ERROR: [$params.atac.assay] The 'qualityThreshold' parameter must be an integer.")
                System.exit(1)
            // Check if the quality threshold is within an acceptable range for a quality score(0-60 inclusive)
            } else if (params.atac.qualityThreshold < 0 || params.atac.qualityThreshold > 60) {
                log.error("ERROR: [$params.atac.assay] The 'qualityThreshold' parameter is not within a valid range.")
                System.exit(1)
            } else {
                value = params.atac.qualityThreshold
            }
        } else {
            value = params.preset.atac.qualityThreshold
        }
        return value
    }

    // Function to check the 'fasta' file parameter
    static <Type> Type validateMergeMethod(params, log) {
        // Initialize map
        def mergeMethodMap = [:]

        // Is the run barcodedTn5?
        if (params.atac.barcodedTn5) {
            Set tinames = params.atac.ti.keySet()

            // Iterate over sample IDs
            for (int j = 0; j < params.atac.sampleIds.size(); j++) {
                String key = params.atac.sampleIds[j]

                // Iterate over tagmentation indexes
                for (int k = 0; k < tinames.size(); k++) {
                    String tikey = tinames[k]
                    Type value

                    if (params.atac.overrides?.get(key)?.get(tikey)?.mergeMethod) {
                        value = params.atac.overrides.get(key).get(tikey).mergeMethod
                    } else if (params.atac.overrides?.get(key)?.mergeMethod) {
                        value = params.atac.overrides.get(key).mergeMethod
                    } else if (params.atac.mergeMethod) {
                        value = params.atac.mergeMethod
                    } else {
                        value = params.preset.atac.mergeMethod
                    }

                    // Performing type checks
                    if (!value.stream().allMatch { it.matches('(both|r1|r2)') }) {
                        log.error("ERROR: [$params.atac.assay] The 'mergeMethod' parameter must be \
                            one of 'both', 'r1', or 'r2'.")
                    } else {
                        mergeMethodMap.put(key + '-' + params.atac.ti.get(tikey), value)
                }
            }
        }
        } else {
            // Iterate over sample IDs
            for (int j = 0; j < params.atac.sampleIds.size(); j++) {
                String key = params.atac.sampleIds[j]
                Type value

                if (params.atac.overrides?.get(key)?.mergeMethod) {
                    value = params.atac.overrides.get(key).mergeMethod
                } else if (params.atac.mergeMethod) {
                    value = params.atac.mergeMethod
                } else {
                    value = params.preset.atac.mergeMethod
                }

                // Performing type checks
                if (!value.stream().allMatch { it.matches('(both|r1|r2)') }) { // groovylint-disable-line
                    log.error("ERROR: [$params.atac.assay] The 'mergeMethod' parameter \
                        must be one of 'both', 'r1', or 'r2'.")
                } else {
                    mergeMethodMap.put(key, value)
            }
    }
        }

        return mergeMethodMap
    }

    // Function for validating crosstalkthreshold and keeping the default value(0.9) if not specified
    static <Type> Type validateCrosstalkthreshold(params, log) {
        // Check if mixed
        if (params.atac.mixed) {
            // Check if set
            if (params.atac.crosstalkthreshold) {
                BigDecimal value = params.atac.crosstalkthreshold
                // Check if set to a value that is valid
                if (1 > value && value > 0.5) {
                    return value
                }
                log.error("ERROR: [$params.atac.assay] The 'crosstalkthreshold' parameter must be a \
                    BigDecimal less than 1 and greater than 0.5.")
                System.exit(1)
            } else {
                BigDecimal value = params.preset.atac.crosstalkthreshold
                return value
            }
        } else {
            return false
        }
    }

    // Function for validating sortSize and keeping the default value(0.25) if not specified
    static <Type> Type validateSortSize(params, log) {
        // Check if set
        if (params.atac.sortSize) {
            BigDecimal value = params.atac.sortSize
            // Check if set to a value that is valid
            if (1 > value && value > 0) {
                return value
            }
            log.error("ERROR: [$params.atac.assay] The 'sortSize' parameter must be a \
                BigDecimal less than 1 and greater than 0.")
            System.exit(1)
        } else {
            BigDecimal value = params.preset.atac.sortSize
            return value
        }
    }

    // Function for validating insert site rounding parameter and keeping default value(0) if not specified
    static <Type> Type validateInsertRounding(params, log) {
        // Check if set
        if (params.atac.rounding) {
            int value = params.atac.rounding
            // Check if set to a value that is valid
            if (Arrays.asList(0, 10, 100, 1000).contains(value)) {
                return value
            }
            log.error("ERROR: [$params.atac.assay] The 'rounding' parameter must be one of 0, 10, 100, or 1000.")
            System.exit(1)
        } else {
            int value = params.preset.atac.rounding
            return value
        }
    }

    // Function for validating max insert size parameter and keeping default value(2000) if not specified
    static <Type> Type validateMaxInsertSize(params, log) {
        // Check if set
        if (params.atac.maxInsertSize) {
            int value = params.atac.maxInsertSize
            // Check if set to a value that is valid, currently the lowest we accept is 100,
            // we could lower this but low values could cause errors in the pipeline
            if (value >= 100) {
                return value
            }
            log.error("ERROR: [$params.atac.assay] The 'maxInsertSize' parameter must be at least 100.")
            System.exit(1)
        } else {
            int value = params.preset.atac.maxInsertSize
            return value
        }
    }

}
