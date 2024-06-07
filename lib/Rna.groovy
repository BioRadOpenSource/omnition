/* groovylint-disable */
import Core
/* groovylint-enable */

class Rna {

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
                    value = params.preset.rna.bead.get(setting)
                }
                // Performing type checks
                if (setting == 'format' && value != null && value.class != String) {
                    log.error("ERROR: [$assayParams.assay] Bead \
                        $setting parameter must be provided as a string: $value")
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

    // Function for setting cell settings and overriding if needed
    static <Type> Type validateCell(params, assayParams, log) {
        // Initializing variables
        List cellSettings = [ 'loaded' ]
        def cellMap = [:]

        // Iterating over setting types
        for (int i = 0; i < cellSettings.size(); i++) {
            String setting = cellSettings[i]
            def sampleMap = [:]

            // Iterating over sample IDs
            for (int j = 0; j < assayParams.sampleIds.size(); j++) {
                String key = assayParams.sampleIds[j]
                Type value

                if (assayParams.overrides?.get(key)?.cell?.get(setting) != null) {
                    value = assayParams.overrides.get(key).cell.get(setting)
                } else if (assayParams.cell?.get(setting) != null) {
                    value = assayParams.cell.get(setting)
                } else {
                    value = params.cell.get(setting)
                }

                // Performing type checks
                if (value != null && value.class != Integer) {
                    log.error("ERROR: [$assayParams.assay] Cell $setting parameter must \
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

    // function to check the 'readQualityScore' parameter
    static Integer validateReadQualityScore(params, assayParams, log) {
        if (params.readQualityScore != null && params.readQualityScore.class != Integer) {
            log.error("ERROR: [$assayParams.assay] The 'quality score'" +
                " parameter must be provided as an integer: $params.readQualityScore")
            System.exit(1)
        } else {
            return assayParams.readQualityScore
        }
    }

    // Function for setting barcode settings and overriding if needed
    static <Type> Type validateBarcode(params, assayParams, log) {
        // Initializing variables
        List barcodeSettings = [ 'force' ]
        def barcodeMap = [:]

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
                    log.error("ERROR: [$assayParams.assay] Barcode $setting \
                        parameter must be provided as an integer: $value")
                    System.exit(1)
                }

                // Add <key, value> pairs to map
                sampleMap.put(key, value)
            }

            // Add <key, value> pairs to map
            barcodeMap.put(setting, sampleMap)
        }

        return barcodeMap
    }

    // Function for validating crosstalkThreshold and keeping the default value(0.9) if not specified
    static <Type> Type validateCrosstalkthreshold(params, assayParams, log) {
        // Check if mixed
        if (assayParams.mixed != null) {
            BigDecimal value
            // Check if set
            if (assayParams.crosstalkThreshold != null) {
                value = assayParams.crosstalkThreshold
                // Check if set to a value that is valid
                if (1 > value && value > 0.5) {
                    return value
                }
                log.error("ERROR: [$assayParams.assay] The 'crosstalkThreshold' \
                    parameter must be a BigDecimal less than 1 and greater than 0.5.")
                System.exit(1)
            } else {
                value = params.preset.rna.crosstalkThreshold
                return (value)
            }
        } else {
            return false
        }
    }

    // Function to check deduplicated bin parameter
    static String validateDedupBin(params, assayParams, log) {
        String value

        // Check if set
        if (assayParams.dedupBin != null) {
            // Check if class type is correct
            if (assayParams.dedupBin.class != Integer) {
                log.error("ERROR: [$assayParams.assay] The 'dedupBin' parameter must be an Integer.")
                System.exit(1)
            } else {
                value = assayParams.dedupBin
            }
        } else {
            value = params.preset.rna.dedupBin
        }
        return value
    }

    // Function for setting the maximum hamming distance for correcting DO UMIs
    static <Type> Type validateUmiHamming(params, assayParams, log) {
        // Initializing variables
        def umiHammingMap = [:]

        // Iterating over sample IDs
        for (int j = 0; j < assayParams.sampleIds.size(); j++) {
            String key = assayParams.sampleIds[j]
            Type value

            if (assayParams.overrides?.get(key)?.umiHamming != null) {
                value = assayParams.overrides.get(key).umiHamming
            } else if (assayParams.umiHamming != null) {
                value = assayParams.umiHamming
            } else {
                value = params.preset.rna.umiHamming
            }

            // Performing type checks
            if (value != null && value.class != Integer) {
                log.error("ERROR: [$assayParams.assay] The 'umiHamming'" +
                    " parameter must be provided as an integer: $value")
                System.exit(1)
            }

            // Add <key, value> pairs to map
            umiHammingMap.put(key, value)
        }

        return umiHammingMap
    }

    /* Function for setting the number of UMIs for considering connections 'real'
        as part of bead merging */
    static <Type> Type validateBeadMergeUmiThreshold(params, assayParams, log) {
        // Initializing variables
        def beadMergeUmiMap = [:]

        // Iterating over sample IDs
        for (int j = 0; j < assayParams.sampleIds.size(); j++) {
            String key = assayParams.sampleIds[j]
            Type value

            if (assayParams.overrides?.get(key)?.beadMergeUmiThreshold != null) {
                value = assayParams.overrides.get(key).beadMergeUmiThreshold
            } else if (assayParams.beadMergeUmiThreshold != null) {
                value = assayParams.beadMergeUmiThreshold
            } else {
                value = params.preset.rna.beadMergeUmiThreshold
            }

            // Performing type checks
            if (value != null && value.class != Integer) {
                log.error("ERROR: [$assayParams.assay] The 'beadMergeUmiThreshold'" +
                    " parameter must be provided as an integer: $value")
                System.exit(1)
            }

            // Add <key, value> pairs to map
            beadMergeUmiMap.put(key, value)
        }

        return beadMergeUmiMap
    }

    // Function to check includeIntrons parameter
    static boolean validateIncludeIntrons(params, assayParams, log) {
        boolean value

        // Check if set
        if (assayParams.includeIntrons != null) {
            // Check if class type is correct
            if (assayParams.includeIntrons.class != Boolean) {
                log.error("ERROR: [$assayParams.assay] The 'includeIntrons' parameter must be boolean (true or false).")
                System.exit(1)
            } else {
                value = assayParams.includeIntrons
            }
        } else {
            value = params.preset.rna.includeIntrons
        }
        return value
    }

    // Function to check includeIntrons parameter
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

}
