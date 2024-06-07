/* groovylint-disable */
import Core
import Atac
import Rna

/* groovylint-enable */

class Validate {

    private final Object workflow
    private final boolean isAWSBatch
    private final Map params
    private final Object runParams
    private final Object log
    private final List messages
    private String assayMode

    Validate(workflow, params, runParams, log, messages) {
        this.workflow = workflow
        this.isAWSBatch = workflow.profile =~ /(awsbatch|tower)/
        this.params = params
        this.runParams = runParams
        this.log = log
        this.messages = messages
        this.setAssayMode()
    }

    //Validations that are always run
    void run() {
        // Limit Tower runs to only the analysis workflows
        // groovylint-disable-next-line
        if ((workflow.profile =~ /(tower)/) && !(workflow.profile =~ /(demo)/) \
        && (runParams.workflow in [ 'full', 'reference' ])) {
            log.error("ERROR: [$runParams.assay] Cannot run the full or reference " \
                        + "workflow on Tower. The 'analysis' workflow is the only " \
                        + "workflow type compatible with Tower. Exiting.")
            exit(1)
        }

        // Validate the output directory
        if (isAWSBatch) {
            Core.validateOutputS3(params, runParams, workflow.profile, log)
        }

        runParams.paramsFile = Core.validateParamsFile(params, assayMode)

        runParams.output = Core.validateOutput(params, runParams, log)
        runParams.resultsDir = "${runParams.output}/Sample_Files/"
        runParams.reportsDir = "${runParams.output}/report/"

        runParams.workflow = Core.validateWorkflow(runParams, log)
        runParams.mixed = Core.validateMixed(runParams, log, messages)
        runParams.reference.directory   = Core.validateReferenceDirectory(runParams, log)

        if (!isAWSBatch) {
            runParams.reference.fasta = Core.validateReferenceFasta(runParams, log)
            runParams.reference.gtf = Core.validateReferenceGtf(runParams, log)
        } else {
            runParams.reference.fasta        = Core.paramFormat(runParams.reference?.fasta)
            runParams.reference.gtf          = Core.paramFormat(runParams.reference?.gtf)
        }
        runParams.species   = Core.getSpecies(runParams.assay, runParams.reference.gtf, log)

        if (assayMode == 'rna') {
            validateRNA()
        }

        if (assayMode == 'atac' || assayMode == 'catac') {
            validateATACAndCATAC()
        }
    }

    //Validations occur only for analysis workflow
    void runAnalysisValidation() {
        // Set workflow params
        if (!isAWSBatch) {
            runParams.input     = Core.validateInput(runParams, log, messages)
        } else {
            runParams.input = Core.validateInputAWS(runParams, runParams.fastqFiles, log, messages)
        }

        runParams.sampleIds = Core.getSampleIds(runParams.assay, runParams.fastqFiles,
             Core.fastqRegEx(), log, messages)

        switch (assayMode) {
            case 'catac' :
                validateCAtacAnalysis()
                break
        }

        if (assayMode == 'atac' || assayMode == 'catac') {
            validateATACAndCATACAnalysis()
        }

        if (assayMode == 'rna') {
            validateRNAAnalysis()
        }
    }

    private void setAssayMode() {
        if (params.atac) {
            assayMode = 'atac'
        }
        if (params.catac) {
            assayMode = 'catac'
        }
        if (params.rna) {
            assayMode = 'rna'
        }
    }

    private void validateRNA() {
        // Ensure that the fasta and gtf are matched
        // Note that the channel can have one or two fasta files in it
        if (!isAWSBatch) { // groovylint-disable-line
            Core.validateReferenceNames(runParams, log)
        }
                // groovylint-disable-next-line
    }

    private void validateATACAndCATAC() {
        runParams.tssWindowSize             = Atac.validateTSSWindow(params, runParams, log)
        runParams.reference.blocklist   = Atac.validateReferenceBlocklist(runParams, log, isAWSBatch)

        // Ensure that the fasta and gtf are matched
        // Note that the channel can have one or two fasta files in it
        Core.validateReferenceNames(runParams, log)
    }

    private void validateRNAAnalysis() {
        runParams.bead               = Rna.validateBead(params, runParams, log)
        runParams.prefix             = Core.validatePrefix(params, runParams, log)
        runParams.barcode            = Rna.validateBarcode(params, runParams, log)
        runParams.crosstalkThreshold = Rna.validateCrosstalkthreshold(params, runParams, log)
        runParams.sampleOverride = Core.validateSampleOverride(runParams, log)
        Rna.validateNumberOfSamples(runParams, messages)
        runParams.beadMergeUmiThreshold = Rna.validateBeadMergeUmiThreshold(params, runParams, log)
        runParams.umiHamming = Rna.validateUmiHamming(params, runParams, log)
        runParams.includeIntrons = Rna.validateIncludeIntrons(params, runParams, log)
        runParams.readQualityScore = Rna.validateReadQualityScore(params, runParams, log)
    }

    private void validateATACAndCATACAnalysis() {
        runParams.prefix              = Core.validatePrefix(params, runParams, log)
        runParams.barcode             = Atac.validateBarcode(params, runParams, log)
        runParams.trimFivePrime       = Atac.validateTrimFivePrime(params, runParams, log)
        runParams.trimThreePrime      = Atac.validateTrimThreePrime(params, runParams, log)
        runParams.sortSize            = Atac.validateSortSize(params, runParams, log)
        runParams.mitoContig          = Atac.validateMitoContig(params, runParams, log)
        runParams.mapQualityThreshold = Atac.validateMapQualityThreshold(params, runParams, log)
        runParams.mergeMethod         = Atac.validateMergeMethod(params, runParams, log)
        runParams.crosstalkThreshold  = Atac.validateCrosstalkthreshold(params, runParams, log)
        runParams.rounding            = Atac.validateInsertRounding(params, runParams, log)
        runParams.maxInsertSize       = Atac.validateMaxInsertSize(params, runParams, log)
        runParams.sampleOverride      = Core.validateSampleOverride(runParams, log)
    }

    private void validateCAtacAnalysis() {
        runParams.ti                 = Atac.validateTI(params, runParams, log, messages)
        runParams.i7AsTi             = Atac.validatei7AsTi(runParams, log, messages)
        runParams.tiRead             = Atac.validateTIRead(params, runParams, log, messages)
        if (isAWSBatch) { // groovylint-disable-line
            if (runParams?.barcodedTn5Config) {
                runParams.barcodedTn5Config = Core.paramFormat(runParams?.barcodedTn5Config) // groovylint-disable-line
            } else {
                runParams.barcodedTn5Config  = Atac.validateBarcodedTn5Config(runParams, log, messages)
            }
        }
        // groovylint-disable-next-line
        runParams.tiErrorOverride    = Atac.validateTIErrorOverride(params, runParams, log, messages)
        runParams.tioverride         = Atac.validateTIOverride(params, runParams, log)
        // groovylint-disable-next-line
    }

}
