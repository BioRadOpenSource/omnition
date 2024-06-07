/* groovylint-disable */
import Core
/* groovylint-enable */

class Atac {

    // Function for setting bead settings and overriding if needed
    static <Type> Type validateBead(params, assayParams, log) {
        // Initializing variables
        List beadSettings = [ 'format']
        def beadMap = [:]

        // Iterating over setting types
        for (int i = 0; i < beadSettings.size(); i++) {
            String setting = beadSettings[i]
            def sampleMap = [:]

            // Iterating over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]
                Type value

                if (assayParams.overrides?.get(key)?.bead?.get(setting) != null) {
                    value = assayParams.overrides.get(key).bead.get(setting)
                } else if (assayParams.bead?.get(setting) != null) {
                    value = assayParams.bead.get(setting)
                } else {
                    value = params.preset.atac.bead.get(setting)
                }

                // Performing type checks
                if (setting == 'format' && value != null && value.class != String) {
                    log.error("ERROR: [$assayParams.assay] Bead $setting parameter \
                        must be provided as a string: $value")
                    System.exit(1)
                }

                // Add <key, value> pairs to map
                sampleMap.put(key, value)
            }

            // Add <key, value> pairs to map
            beadMap.put(setting, sampleMap)
        }

        return beadMap
    }

    // Function to check the 'tiErrorOverride' parameter
    static boolean validateTIErrorOverride(params, assayParams, log, messages) {
        boolean value

        // Check if set
        if (assayParams.tiErrorOverride != null) {
            // Check if class type is correct
            if (assayParams.tiErrorOverride.class != Boolean) {
                log.error("ERROR: [$assayParams.assay] The 'tiErrorOverride' parameter must be boolean(true or false).")
                System.exit(1)
            } else {
                log.warn("[$assayParams.assay] Overriding errors in the config file. Only TIs specified \
                    in the config file will be used.")
                messages.add("WARN: [$assayParams.assay] Overriding errors in the config file. \
                    Only TIs specified in the config file will be used.")
                value = assayParams.tiErrorOverride
            }
        } else {
            value = params.preset.atac.tiErrorOverride
        }
        return value
    }

    // Function to check the barcodedTn5 config parameter
    static <Type> Type validateBarcodedTn5Config(params, log, messages) {
        Type value

        // Check if cATAC assay
        if (params.assay != null && params.assay.toLowerCase() == "catac") {
            // Check if parameter set
            if (params.barcodedTn5Config != null) {
                value = params.barcodedTn5Config
                File file = new File(value)
                // Check if file exists
                if (file.exists()) {
                    // Check if file is a csv
                    if (value.stream().allMatch { it.matches('\\S+(.csv)') }) {
                        validateBarcodedTn5ConfigFile(file, log)
                        return value
                    }
                    log.error("ERROR: [$params.assay] The 'BarcodedTn5Config' parameter file \
                        must be a csv file $value")
                    System.exit(1)
                } else {
                    log.error("ERROR: [$params.assay] The 'BarcodedTn5Config' parameter \
                        file does not exist: $value")
                    System.exit(1)
                }
            } else {
                log.info("INFO: [$params.assay] No 'BarcodedTn5Config' Config file provided, \
                    executing superloading workflow.")
                messages.add("INFO: [$params.assay] No 'BarcodedTn5Config' Config file provided, \
                    executing superloading workflow.")
                return false
            }
        } else {
            return false
        }
    }

    static void validateBarcodedTn5ConfigFile(file, log) {
        List lines = file.readLines()
        for (int i = 1; i < lines.size(); i++) {
            List fields = lines[i].split(',')
            if (!fields[0].matches('[a-zA-Z0-9_-]+')) {
                log.error("ERROR: [$file] The file contains a sample_name with an invalid character")
                System.exit(1)
            }
        }
    }

    // Function to check the TI parameter
    static <Type> Type validateTI(params, assayParams, log, messages) {
        Type value

        // Check if catac workflow
        if (assayParams.assay != null && assayParams.assay.toLowerCase() == "catac") { // groovylint-disable-line
            // Check if parameter set
            if (assayParams.ti != null) {
                value = assayParams.ti
            } else {
                log.info("INFO: [$assayParams.assay] No TI list provided, using preset TIs.")
                messages.add("INFO: [$assayParams.assay] No TI list provided, using preset TIs.")
                value = params.preset.atac.ti
            }

            return(value)
        }
        return false
    }

    // Function to check the TIread parameter
    static <Type> Type validateTIRead(params, assayParams, log, messages) {
        Type value

        // Check if catac workflow
        if (assayParams.assay != null && assayParams.assay.toLowerCase() == "catac") { // groovylint-disable-line
            // Check if parameter set
            if (assayParams.tiRead != null) {
                value = assayParams.tiRead

                // Check if i7AsTi is false
                if (assayParams.i7AsTi == true) {
                    log.error("ERROR: [$assayParams.assay] The 'tiRead' parameter cannot be used in \
                        conjuction with the 'i7AsTi' parameter.")
                }

                // Performing type checks
                if (!value.matches('(r1|r2)')) {
                    log.error("ERROR: [$assayParams.assay] The 'tiRead' parameter must be 'r1', or 'r2'.")
                }
            } else {
                value = params.preset.atac.tiRead
            }

            if (assayParams.i7AsTi == true) {
                log.info("INFO: [$assayParams.assay] TIs are set to be located in the i7 barcode.")
                messages.add("INFO: [$assayParams.assay] TIs are set to be located in the i7 barcode.")
            } else {
                log.info("INFO: [$assayParams.assay] TIs are set to be located in $value.")
                messages.add("INFO: [$assayParams.assay] TIs are set to be located in $value.")
            }
            return(value)
        }
        return false
    }

    // Function to check the 'i7AsTi' parameter
    static boolean validatei7AsTi(params, log, messages) {
        boolean value

        // Check if set
        if (params.i7AsTi != null) {
            // Check if class type is correct
            if (params.i7AsTi.class != Boolean) {
                log.error("ERROR: [$params.assay] The 'i7AsTi' parameter must be boolean(true or false).")
                System.exit(1)
            } else if (params.assay != "cATAC") { // groovylint-disable-line
                log.error("ERROR: [$params.assay] The assay must be the 'cATAC' assay \
                    in order to set 'i7AsTi' parameter to true.")
                System.exit(1)
            } else {
                String message = "INFO: [$params.assay] Using i7 sequence as a tagmentation index."
                log.info(message)
                messages.add(message)
                value = params.i7AsTi
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
        if (params.config != null) {
            value = params.config
            // Check if file exists
            if (new File(value).exists()) {
                return value
            }
            log.error("ERROR: [$params.assay] The 'config' parameter file does not exist: $value")
            System.exit(1)
        } else {
            Core.paramError(params.assay, 'config', log)
        }
    }

    // Function for setting the overrides for specific TIs if needed
    static boolean validateTIOverride(assayParams, params, log) {
        // Initializing variables
        boolean value
        // Iterating over the provided override sample IDs
        if (params.overrides != null && assayParams.assay != null && assayParams.assay.toLowerCase() == "catac") { // groovylint-disable-line
            Set overridesamples = params.overrides.keySet()

            for (int i = 0; i < overridesamples.size(); i++) {
                String key = overridesamples[i]
                Set tinames = params.ti.keySet()
                Set fastqLevelOverrides = ["barcode", "downsample"]
                // Get parameter overrides for the sample
                def entries = params.overrides.get(key)
                .collectEntries { entry -> [ entry.key, entry.value ] }
                for (e in entries) {
                    // If the override is a TI, class will be a Map not a string, where the key is the TI
                    if (e.value.getClass() != String) {
                        if (tinames.contains(e.key) || fastqLevelOverrides.contains(e.key)) { // Valid TI
                            value = true
                        } else {
                            log.error("ERROR: [$params.assay] The 'overrides' parameter \
                                $e.key is neither in the ti list: $tinames nor is it a recognized \
                                fastq level override.")
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

    // Function to check reference Blocklist
    static <Type> Type validateReferenceBlocklist(params, log, isAWSBatch) {
        Type value
        int count

        // Check if set
        if (params.reference?.blocklist != null) {
            value = Core.paramFormat(params.reference?.blocklist)
            // Check if endings are valid
            if (value.stream().allMatch { it.matches('\\S+(blocklist.bed)') }) {
                //If we are using AWS Batch mode, grab the files from the blocklistFiles
                if (isAWSBatch) {
                    value = params.reference.blocklistFiles
                }
                if (value.stream().allMatch {
                    (!isAWSBatch && new File(it).exists()) || (isAWSBatch && it.exists()) }) {
                    count = Core.countElement(value)
                    // Check if valid number of files have been provided
                    if (params.mixed != true && count != 1) {
                        log.error("ERROR: [$params.assay] The reference 'blocklist' parameter accepts at most \
                            one file per reference genome. If using two references, set \
                            'mixed: true' in the parameters file.")
                        System.exit(1)
                    } else if (params.mixed == true && count > 2) {
                        log.error("ERROR: [$params.assay] The reference 'blocklist' parameter accepts \
                            at most one file per reference genome.")
                        System.exit(1)
                    } else {
                        return value
                    }
                } else {
                    log.error("ERROR: [$params.assay] The reference 'blocklist' parameter \
                        file does not exist: $value")
                    System.exit(1)
                }
            } else {
                log.error("ERROR: [$params.assay] The reference 'blocklist' \
                    parameters must end in 'blocklist.bed'.")
                System.exit(1)
            }
        } else {
            // Check if workflow is analysis
            if (params.workflow == 'analysis') {
                Core.paramError(params.assay + '_' + params.workflow, 'blocklist', log)
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
            log.error("ERROR: [$assay] No FASTQ files found in input directory. \
                Check parameters and file name requirements.")
            System.exit(1)
        } else {
            return sampleIdList
        }
    }

    // Function for setting barcode settings and overriding if needed
    static <Type> Type validateBarcode(params, assayParams, log) {
        // Initializing variables
        List barcodeSettings = [ 'force' ]
        def barcodeMap = [:]

        // Is the run cATAC?
        if (assayParams.assay != null && assayParams.assay.toLowerCase() == "catac") { // groovylint-disable-line
            Set tinames = assayParams.ti.keySet()

            // Iterating over setting types
            for (int i = 0; i < barcodeSettings.size(); i++) {
                String setting = barcodeSettings[i]
                def sampleMap = [:]

                // Iterating over sample IDs
                for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                    String key = assayParams.sampleIds[j]

                    // Iterating over tagmentation indexes
                    for (int k = 0; k < tinames.size(); k++) {
                        String tikey = tinames[k]
                        Type value

                        if (assayParams.overrides?.get(key)?.get(tikey)?.barcode?.get(setting) != null) {
                            value = assayParams.overrides.get(key).get(tikey).barcode.get(setting)
                        } else if (assayParams.overrides?.get(key)?.barcode?.get(setting) != null) {
                            value = assayParams.overrides.get(key).barcode.get(setting)
                        } else if (assayParams.barcode?.get(setting) != null) {
                            value = assayParams.barcode.get(setting)
                        } else {
                            value = params.barcode.get(setting)
                        }

                        // Performing type checks
                        if (value != null && value.class != Integer) {
                            log.error("ERROR: [$assayParams.assay] Barcode $setting parameter must \
                                be provided as an integer: $value")
                            System.exit(1)
                        }

                        // Add <key, value> pairs to map
                        sampleMap.put(key + '-' + assayParams.ti.get(tikey), value)
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
                for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                    String key = assayParams.sampleIds[j]
                    Type value

                    if (assayParams.overrides?.get(key)?.barcode?.get(setting) != null) {
                        value = assayParams.overrides.get(key).barcode.get(setting)
                    } else if (assayParams.barcode?.get(setting) != null) {
                        value = assayParams.barcode.get(setting)
                    } else {
                        value = params.barcode.get(setting)
                    }

                    // Performing type checks
                    if (value != null && value.class != Integer) {
                        log.error("ERROR: [$assayParams.assay] Barcode $setting parameter must be \
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

    // Function for setting the number of bases to trim from 5' end of R2 reads
    static <Type> Type validateTrimFivePrime(params, assayParams, log) {
        // Initializing variables
        def trimMap = [:]

        // Is the run cATAC?
        if (assayParams.assay != null && assayParams.assay.toLowerCase() == "catac") { // groovylint-disable-line
            Set tinames = assayParams.ti.keySet()

            // Iterating over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]

                // Iterating over tagmentation indexes
                for (int k = 0; k < tinames.size(); k++) {
                    String tikey = tinames[k]
                    Type value
                    if (assayParams.overrides?.get(key)?.get(tikey)?.trimFivePrime != null) {
                        value = assayParams.overrides.get(key).get(tikey).trimFivePrime
                    } else if (assayParams.overrides?.get(key)?.trimFivePrime != null) {
                        value = assayParams.overrides.get(key).trimFivePrime
                    } else if (assayParams.trimFivePrime != null) {
                        value = assayParams.trimFivePrime
                    } else {
                        value = params.preset.atac.trimFivePrime
                    }

                    // Performing type checks
                    if (value != null && value.class != Integer) {
                        log.error("ERROR: [$assayParams.assay] The 'trimFivePrime'" +
                            " parameter must be provided as an integer: $value")
                        System.exit(1)
                    }

                    // Add <key, value> pairs to map
                    trimMap.put(key + '-' + assayParams.ti.get(tikey), value)
                }
            }
        } else {
            // Iterating over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]
                Type value

                if (assayParams.overrides?.get(key)?.trimFivePrime != null) {
                    value = assayParams.overrides.get(key).trimFivePrime
                } else if (assayParams.trimFivePrime != null) {
                    value = assayParams.trimFivePrime
                } else {
                    value = params.preset.atac.trimFivePrime
                }

                // Performing type checks
                if (value != null && value.class != Integer) {
                    log.error("ERROR: [$assayParams.assay] The 'trimFivePrime'" +
                        " parameter must be provided as an integer: $value")
                    System.exit(1)
                }

                // Add <key, value> pairs to map
                trimMap.put(key, value)
            }
        }

        return trimMap
    }

    // Function for setting the number of bases to trim from 3' end of R2 reads
    static <Type> Type validateTrimThreePrime(params, assayParams, log) {
        // Initializing variables
        def trimMap = [:]

        // Is the run cATAC?
        if (assayParams.assay != null && assayParams.assay.toLowerCase() == "catac") { // groovylint-disable-line
            Set tinames = assayParams.ti.keySet()

            // Iterating over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]

                // Iterating over tagmentation indexes
                for (int k = 0; k < tinames.size(); k++) {
                    String tikey = tinames[k]
                    Type value
                    if (assayParams.overrides?.get(key)?.get(tikey)?.trimThreePrime != null) {
                        value = assayParams.overrides.get(key).get(tikey).trimThreePrime
                    } else if (assayParams.overrides?.get(key)?.trimThreePrime != null) {
                        value = assayParams.overrides.get(key).trimThreePrime
                    } else if (assayParams.trimThreePrime != null) {
                        value = assayParams.trimThreePrime
                    } else {
                        value = params.preset.atac.trimThreePrime
                    }

                    // Performing type checks
                    if (value != null && value.class != Integer) {
                        log.error("ERROR: [$assayParams.assay] The 'trimThreePrime'" +
                            " parameter must be provided as an integer: $value")
                        System.exit(1)
                    }

                    // Add <key, value> pairs to map
                    trimMap.put(key + '-' + assayParams.ti.get(tikey), value)
                }
            }
        } else {
            // Iterating over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]
                Type value

                if (assayParams.overrides?.get(key)?.trimThreePrime != null) {
                    value = assayParams.overrides.get(key).trimThreePrime
                } else if (assayParams.trimThreePrime != null) {
                    value = assayParams.trimThreePrime
                } else {
                    value = params.preset.atac.trimThreePrime
                }

                // Performing type checks
                if (value != null && value.class != Integer) {
                    log.error("ERROR: [$assayParams.assay] The 'trimThreePrime'" +
                        " parameter must be provided as an integer: $value")
                    System.exit(1)
                }

                // Add <key, value> pairs to map
                trimMap.put(key, value)
            }
        }

        return trimMap
    }

    // Function to check the 'mitoContig' parameter
    static <Type> Type validateMitoContig(params, assayParams, log) {
        String value

        // Check if set
        if (assayParams.mitoContig != null) {
            // Check if class type is correct
            if (assayParams.mitoContig.class != String) {
                log.error("ERROR: [$assayParams.assay] The 'mitoContig' parameter must be a string.")
                System.exit(1)
            } else {
                value = assayParams.mitoContig
            }
        } else {
            value = params.preset.atac.mitoContig
        }
        return value
    }

    // Function to check the 'tssWindowSize' parameter
    static <Type> Type validateTSSWindow(params, assayParams, log) {
        int value

        // Check if set
        if (assayParams.tssWindowSize != null) {
            // Check if class type is correct
            if (assayParams.tssWindowSize.class != Integer) {
                log.error("ERROR: [$assayParams.assay] The 'tssWindowSize' parameter must be an even integer.")
                System.exit(1)
            } else {
                // Check if given value is even
                if (assayParams.tssWindowSize % 2 == 0) {
                    value = assayParams.tssWindowSize
                } else {
                    log.error("ERROR: [$assayParams.assay] The 'tssWindowSize' parameter must be an even integer.")
                }
            }
        } else {
            value = params.preset.atac.tssWindowSize
        }
        return value
    }

    // Function to check the 'mapQualityThreshold' parameter
    static <Type> Type validateMapQualityThreshold(params, assayParams, log) {
        int value
        // Check if set
        if (assayParams.mapQualityThreshold != null) {
            // Check if class type is correct
            if (assayParams.mapQualityThreshold.class != Integer) {
                log.error("ERROR: [$assayParams.assay] The 'mapQualityThreshold' parameter must be an integer.")
                System.exit(1)
            // Check if the quality threshold is within an acceptable range for a quality score(0-60 inclusive)
            } else if (assayParams.mapQualityThreshold < 0 || assayParams.mapQualityThreshold > 60) {
                log.error("ERROR: [$assayParams.assay] The 'mapQualityThreshold' \
                    parameter is not within a valid range.")
                System.exit(1)
            } else {
                value = assayParams.mapQualityThreshold
            }
        } else {
            value = params.preset.atac.mapQualityThreshold
        }
        return value
    }

    // Function to check the 'fasta' file parameter
    static <Type> Type validateMergeMethod(params, assayParams, log) {
        // Initialize map
        def mergeMethodMap = [:]

        // Is the run cATAC?
        if (assayParams.assay != null && assayParams.assay.toLowerCase() == "catac") { // groovylint-disable-line
            Set tinames = assayParams.ti.keySet()

            // Iterate over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]

                // Iterate over tagmentation indexes
                for (int k = 0; k < tinames.size(); k++) {
                    String tikey = tinames[k]
                    Type value

                    if (assayParams.overrides?.get(key)?.get(tikey)?.mergeMethod != null) {
                        value = assayParams.overrides.get(key).get(tikey).mergeMethod
                    } else if (assayParams.overrides?.get(key)?.mergeMethod != null) {
                        value = assayParams.overrides.get(key).mergeMethod
                    } else if (assayParams.mergeMethod != null) {
                        value = assayParams.mergeMethod
                    } else {
                        value = params.preset.atac.mergeMethod
                    }

                    // Performing type checks
                    if (!value.stream().allMatch { it.matches('(both|r1|r2)') }) {
                        log.error("ERROR: [$assayParams.assay] The 'mergeMethod' parameter must be \
                            one of 'both', 'r1', or 'r2'.")
                    } else {
                        mergeMethodMap.put(key + '-' + assayParams.ti.get(tikey), value)
                    }
                }
            }
        } else {
            // Iterate over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]
                Type value

                if (assayParams.overrides?.get(key)?.mergeMethod != null) {
                    value = assayParams.overrides.get(key).mergeMethod
                } else if (assayParams.mergeMethod != null) {
                    value = assayParams.mergeMethod
                } else {
                    value = params.preset.atac.mergeMethod
                }

                // Performing type checks
                if (!value.stream().allMatch { it.matches('(both|r1|r2)') }) { // groovylint-disable-line
                    log.error("ERROR: [$assayParams.assay] The 'mergeMethod' parameter \
                        must be one of 'both', 'r1', or 'r2'.")
                } else {
                    mergeMethodMap.put(key, value)
                }
            }
        }

        return mergeMethodMap
    }

    // Function for validating crosstalkThreshold and keeping the default value(0.9) if not specified
    static <Type> Type validateCrosstalkthreshold(params, assayParams, log) {
        // Check if mixed
        if (assayParams.mixed != null) {
            // Check if set
            if (assayParams.crosstalkThreshold != null) {
                BigDecimal value = assayParams.crosstalkThreshold
                // Check if set to a value that is valid
                if (1 > value && value > 0.5) {
                    return value
                }
                log.error("ERROR: [$assayParams.assay] The 'crosstalkThreshold' parameter must be a \
                    BigDecimal less than 1 and greater than 0.5.")
                System.exit(1)
            } else {
                BigDecimal value = params.preset.atac.crosstalkThreshold
                return value
            }
        } else {
            return false
        }
    }

    // Function for validating sortSize and keeping the default value(0.25) if not specified
    static <Type> Type validateSortSize(params, assayParams, log) {
        // Check if set
        if (assayParams.sortSize != null) {
            BigDecimal value = assayParams.sortSize
            // Check if set to a value that is valid
            if (1 > value && value > 0) {
                return value
            }
            log.error("ERROR: [$assayParams.assay] The 'sortSize' parameter must be a \
                BigDecimal less than 1 and greater than 0.")
            System.exit(1)
        } else {
            BigDecimal value = params.preset.atac.sortSize
            return value
        }
    }

    // Function for validating insert site rounding parameter and keeping default value(0) if not specified
    static <Type> Type validateInsertRounding(params, assayParams, log) {
        // Check if set
        if (assayParams.rounding != null) {
            int value = assayParams.rounding
            // Check if set to a value that is valid
            if (Arrays.asList(0, 10, 100, 1000).contains(value)) {
                return value
            }
            log.error("ERROR: [$assayParams.assay] The 'rounding' parameter must be one of 0, 10, 100, or 1000.")
            System.exit(1)
        } else {
            int value = params.preset.atac.rounding
            return value
        }
    }

    // Function for validating max insert size parameter and keeping default value(2000) if not specified
    static <Type> Type validateMaxInsertSize(params, assayParams, log) {
        // Check if set
        if (assayParams.maxInsertSize != null) {
            int value = assayParams.maxInsertSize
            // Check if set to a value that is valid, currently the lowest we accept is 100,
            // we could lower this but low values could cause errors in the pipeline
            if (value >= 100) {
                return value
            }
            log.error("ERROR: [$assayParams.assay] The 'maxInsertSize' parameter must be at least 100.")
            System.exit(1)
        } else {
            int value = params.preset.atac.maxInsertSize
            return value
        }
    }

}
