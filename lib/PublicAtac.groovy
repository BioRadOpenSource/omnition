
class PublicAtac {

    // Function for setting barcode settings and overriding if needed
    static <Type> Type validateBarcodeNfSchema(assayParams) {
        // Initializing variables
        def barcodeMap = [:]
        def sampleMap = [:]

        for (int j = 0; j < assayParams.tiMetaMap.size(); j++) {
            def currentMap = assayParams.tiMetaMap[j][0]
            for (key in currentMap.keySet() - "sampleId") { // Exclude 'sampleId' key
                Type value
                String sampleId = currentMap['sampleId'] // groovylint-disable-line
                if (currentMap.get(key)?.get('barcode') != []) { // groovylint-disable-line
                    value = currentMap.get(key).get('barcode').get('force') // groovylint-disable-line
                } else {
                    // Using the sample-level default value
                    value = assayParams.get('barcode').get('force').get(sampleId) // groovylint-disable-line
                }
                // Add <key, value> pairs to map
                sampleMap.put(sampleId + '-' + currentMap.get(key).get('sequence'), value)
            }
        }
        // Add <key, value> pairs to map
        barcodeMap.put('force', sampleMap) // groovylint-disable-line
        return barcodeMap
    }

    // Function for setting the number of bases to trim from the 5' or 3' end of R2 reads
    static <Type> Type validateTrimFiveOrThreePrimeNfSchema(endToTrrim, assayParams) {
        // Initializing variables
        def trimMap = [:]
        // Iterating over sample IDs
        for (int j = 0; j < assayParams.tiMetaMap.size(); j++) {
            def currentMap = assayParams.tiMetaMap[j][0]
            for (key in currentMap.keySet() - "sampleId") { // groovylint-disable-line
                Type value
                String sampleId = currentMap['sampleId'] // groovylint-disable-line
                if (currentMap.get(key)?.get(endToTrrim) != []) {
                    value = currentMap.get(key).get(endToTrrim)
                } else {
                    // Using the sample-level default value
                    value = assayParams.get(endToTrrim).get(sampleId)
                }
                // Add <key, value> pairs to map
                trimMap.put(sampleId + '-' + currentMap.get(key).get('sequence'), value)
            }
        }
        return trimMap
    }

}
