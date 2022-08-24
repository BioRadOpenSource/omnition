#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

log.info(Core.bioradLogo(params.monochrome_logs))
messages = []

// Check if help flag was given
if (params.help) {
    message = Core.helpMessage(params.monochrome_logs)
    log.info(message)
    System.exit(0)
}

// Validate the output directory
Core.validateOutputDir(params, log)

if (params.atac) {
    // Initialize mapped parameters if they were not set
    params.atac.reference   = params.atac.reference ?: [:]
    params.atac.contaminant = params.atac.contaminant ?: [:]

    // Set global assay params
    params.atac.assay                     = 'ATAC'
    params.atac.workflow                  = Atac.validateWorkflow(params, log)
    params.atac.mixed                     = Atac.validateMixed(params, log, messages)
    params.atac.tssWindowSize             = Atac.validateTSSWindow(params, log)
    params.atac.reference.directory   = Atac.validateReferenceDirectory(params, log)
    params.atac.contaminant.directory = Atac.validateContaminantDirectory(params, log)
    params.atac.reference.fasta       = Atac.validateReferenceFasta(params, log)
    params.atac.reference.gtf         = Atac.validateReferenceGtf(params, log)
    params.atac.reference.blocklist   = Atac.validateReferenceBlocklist(params, log)

    // Ensure that the fasta and gtf are matched
    // Note that the channel can have one or two fasta files in it
    Core.validateReferenceNames(params.atac, log)

    // Validate reference workflow
    if (params.atac.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.atac.assay] Executing reference workflow.")
        messages.add("INFO: [$params.atac.assay] Executing reference workflow.")
    }

    // Validate analysis workflow
    if (params.atac.workflow in [ 'analysis', 'full' ]) {
        log.info("INFO: [$params.atac.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.atac.assay] Executing analysis workflow.")

        // Set workflow params
        params.atac.input              = Atac.validateInput(params, log, messages)
        params.atac.sampleIds          = Atac.getSampleIds(params.atac.assay, params.atac.input,
            [ '*R{1,2}_*.{fq,fastq}.gz', '*_{1,2}.{fq,fastq}.gz'], log) // groovylint-disable-line
        params.atac.species            = Core.getSpecies(params.atac.assay, params.atac.reference.gtf, log)
        params.atac.barcodedTn5        = Atac.validateBarcodedTn5(params, log, messages)
        params.atac.ti                 = Atac.validateTI(params, log, messages)
        params.atac.i7asti             = Atac.validatei7asti(params, log, messages)
        params.atac.tiread             = Atac.validateTIRead(params, log, messages)
        params.atac.cell               = Atac.validateCell(params, log)
        params.atac.barcode            = Atac.validateBarcode(params, log)
        params.atac.trim               = Atac.validateTrim(params, log)
        params.atac.sortSize           = Atac.validateSortSize(params, log)
        params.atac.mitoContig         = Atac.validateMitoContig(params, log)
        params.atac.qualityThreshold   = Atac.validateQualityThreshold(params, log)
        params.atac.mergeMethod        = Atac.validateMergeMethod(params, log)
        params.atac.crosstalkthreshold = Atac.validateCrosstalkthreshold(params, log)
        params.atac.rounding           = Atac.validateInsertRounding(params, log)
        params.atac.maxInsertSize      = Atac.validateMaxInsertSize(params, log)
        params.atac.barcodedTn5Config  = Atac.validateBarcodedTn5Config(params, log, messages)
        params.atac.tierroroverride    = Atac.validateTIErrorOverride(params, log, messages)
        params.atac.sampleoverride     = Atac.validateSampleOverride(params, log)
        params.atac.tioverride         = Atac.validateTIOverride(params, log)
    }
} else {
    log.error('ERROR: No assay(s) specified. Check parameters file.')
    exit(1)
}

// Import workflows
include { PULL_IMAGES } from './workflows/core/singularity.nf'
include { ATAC_ANALYSIS } from './workflows/atac/atac_analysis.nf'
include { ATAC_REFERENCE } from './workflows/atac/atac_reference.nf'

workflow {
    // Analyze ATAC data
    if (params.atac) {
        if (workflow.containerEngine == 'singularity') {
            workflow_modules = Channel.fromPath("${baseDir}/modules/atac/*.nf")
            core_modules = Channel.fromPath("${baseDir}/modules/core/*.nf")

            PULL_IMAGES(
                workflow_modules.mix(core_modules)
           )
        }

        ch_images_pulled = workflow.containerEngine == 'singularity' ?
            PULL_IMAGES.out.images_pulled : Channel.of(true).first()

        // Create a file in the results direcotry of all messages passed during validation
        ch_atac_messages  = Channel.fromList(messages).collectFile(
            name: "${params.reportsDir}/messages.txt", newLine: true, sort: 'index')

        // Generate references
        if (params.atac.workflow in [ 'reference', 'full' ]) {
            // Set workflow channels
            ch_atac_reference_fasta = Channel.fromPath(params.atac.reference.fasta)
            ch_atac_reference_gtf   = Channel.fromPath(params.atac.reference.gtf)

            // Execute reference workflow
            ATAC_REFERENCE(
                ch_atac_reference_fasta,
                ch_atac_reference_gtf,
                ch_images_pulled
           )
        }

        // Analyze data
        if (params.atac.workflow in [ 'analysis', 'full' ]) {
            // Set workflow channels
            ch_atac_fastq  = Channel.fromFilePairs([ "${params.atac.input}/*_{1,2}.{fq,fastq}.gz",
                "${params.atac.input}/*_R{1,2}_*.{fq,fastq}.gz" ], flat: true)
                .map { prefix, file1, file2 -> tuple(prefix - ~/(_L[0-9][0-9][0-9])$/, file1, file2) }
                .groupTuple()
                .filter { !(it =~ /Undetermined/) }
                .ifEmpty {
                    log.error("ERROR: [$params.atac.assay] No FASTQ read files in input directory. \
                    Check parameters and/or see documentation for file naming guidelines.")
                    exit(1)
                }
            ch_atac_TI_config = params.atac.barcodedTn5Config ? Channel.fromPath(params.atac.barcodedTn5Config) :
                Channel.empty()
            if (params.atac.mixed) {
                ch_atac_reference_fasta = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.fasta      :
                    Channel.fromPath("${params.atac.reference.directory}/mixed.fa").first()
                ch_atac_reference_gtf   = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.gtf        :
                    Channel.fromPath("${params.atac.reference.directory}/mixed.filtered.gtf").first()
            } else {
                ch_atac_reference_fasta = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.fasta  :
                    Channel.fromPath(params.atac.reference.fasta).first()
                ch_atac_reference_gtf   = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.gtf    :
                    Channel.fromPath(params.atac.reference.gtf).first()
            }
            ch_atac_reference_blocklist = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.blocklist  :
                Channel.fromPath(params.atac.reference.blocklist).first()
            ch_atac_reference_index     = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.bwaIndex   :
                Channel.fromPath("${params.atac.reference.directory}/bwa-index/").first()
            ch_atac_reference_sizes     = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.genomeSize :
                Channel.fromPath("${params.atac.reference.directory}/genome.sizes").first()
            ch_atac_archr_ref           = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.archrRef   :
                Channel.fromPath("${params.atac.reference.directory}/archr/BSgenome.ref.na.1.0_1.0.tar.gz").first()
            ch_atac_tss                 = params.atac.workflow == 'full' ? ATAC_REFERENCE.out.tss        :
                Channel.fromPath("${params.atac.reference.directory}/TSS.${params.atac.tssWindowSize}.bed").first()

            // Catch missing reference files
            Channel.empty()
                .mix(ch_atac_reference_fasta,
                    ch_atac_reference_gtf,
                    ch_atac_reference_index,
                    ch_atac_reference_sizes,
                    ch_atac_archr_ref,
                    ch_atac_reference_blocklist,
                    ch_atac_tss)
                .flatten()
                .map { item -> [ item, item.exists() ] }
                .filter { !it[1] }
                .collect { it[0] }
                .subscribe {
                    log.error('ERROR: [ATAC] Missing ATAC reference file(s). Check parameters and/or run\
                        reference workflow.\n  ' + it.join('\n  ')) // groovylint-disable-line
                    exit(1)
                }

            // Execute workflow
            ATAC_ANALYSIS(
                ch_atac_fastq,
                ch_atac_reference_fasta,
                ch_atac_reference_gtf,
                ch_atac_reference_index,
                ch_atac_reference_sizes,
                ch_atac_archr_ref,
                ch_atac_reference_blocklist,
                ch_atac_tss,
                ch_atac_TI_config,
                ch_atac_messages,
                params.atac.assay,
                ch_images_pulled
           )
        }
    }
}
