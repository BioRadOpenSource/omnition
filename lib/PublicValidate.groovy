/* groovylint-disable */
import PublicCore
import PublicAtac

/* groovylint-enable */

class PublicValidate {

    private final Object workflow
    private final boolean isAWSBatch
    private final Map params
    private final Object runParams
    private final Map coreParams
    private final Object log
    private final List messages
    private final String assayMode

    /* groovylint-disable */
    PublicValidate(workflow, params, runParams, coreParams, log, messages, assayMode) {
        this.workflow = workflow
        this.isAWSBatch = workflow.profile =~ /(awsbatch|tower)/
        this.params = params
        this.runParams = runParams
        this.coreParams = coreParams
        this.log = log
        this.messages = messages
        this.assayMode = assayMode
    }
    /* groovylint-enable */

    //Validations that are always run
    void run() {
        if ((workflow.profile =~ /(tower)/) && !(workflow.profile =~ /(demo)/) \
        && !(runParams.input =~ /(automated-testing)/) \
        && (runParams.workflow in [ 'full', 'reference' ])) {
            log.error("ERROR: [$runParams.assay] Cannot run the full or reference " \
                        + "workflow on Tower. The 'analysis' workflow is the only " \
                        + "workflow type compatible with Tower. Exiting.")
            exit(1)
        }

        runParams.paramsFile = PublicCore.validateParamsFile(params)

        runParams.outputDir = PublicCore.validateOutput(params, coreParams, log)
        runParams.resultsDir = "${runParams.outputDir}/Sample_Files/"
        runParams.reportsDir = "${runParams.outputDir}/report/"

        runParams.workflow = PublicCore.validateWorkflow(runParams, log)
        runParams.reference.directory = PublicCore.validateReferenceDirectory(runParams, log)
        // Cytometry doesn't use fasta or gtf files so we don't need to validate them
        if (assayMode != "cytometry") { // groovylint-disable-line
            runParams.species = PublicCore.getSpeciesSchema(runParams.assay, runParams.reference.gtf, log)
        } else {
            runParams.species = null
        }
        PublicCore.loadSampleSheetParams(runParams)

        switch (assayMode) {
            case 'atac':
                break
            case 'catac':
                break
            case 'rna':
                break
        }
    }

    // Validations occur only for analysis workflow
    void runAnalysisValidation() {
        // Set workflow params
        if (!isAWSBatch) {
            runParams.input = PublicCore.validateInput(runParams, log, messages)
        } else {
            runParams.input = PublicCore.validateInputAWS(runParams, params.core.fastqFiles, log, messages)
        }

        String fastqRegex = PublicCore.fastqRegEx() // groovylint-disable-line

        runParams.sampleIds          = PublicCore.getSampleIds(runParams.assay, params.core.fastqFiles,
             fastqRegex, log, messages)
        runParams.sampleIds          = PublicCore.updateSampleIds(runParams)
        PublicCore.validateNumberOfSamples(runParams, messages)
        if (!runParams.sampleSheetMetaMap) {
            PublicCore.loadSampleLevelParams(runParams)
        }

        if (assayMode == 'catac') {
            runParams.trimThreePrime = PublicAtac.validateTrimFiveOrThreePrimeNfSchema("trimThreePrime", runParams)
            runParams.trimFivePrime = PublicAtac.validateTrimFiveOrThreePrimeNfSchema("trimFivePrime", runParams)
            runParams.barcode = PublicAtac.validateBarcodeNfSchema(runParams)
        }

        switch (assayMode) {
            case 'catac':
                break
            case 'atac':
                break
            case 'rna':
                break
            case 'cytometry':
                break
            default:
                throw new IllegalArgumentException("Unknown assay mode: " + assayMode)
        }
    }

}
