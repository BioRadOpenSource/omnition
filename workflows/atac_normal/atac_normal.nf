/*
Validate parameters and execute atac normal workflows
*/

// Import workflows
include { PULL_IMAGES } from '../core/singularity.nf'
include { ATAC_NORMAL_ANALYSIS } from './atac_normal_analysis.nf'
include { ATAC_NORMAL_REFERENCE } from './atac_normal_reference.nf'
include { GET_FASTQ_FILES } from '../core/fastq_files_array.nf'
include { GET_BLOCKLIST_FILES } from '../core/blocklist_files.nf'

workflow ATAC_NORMAL {
    take:
    messages

    main:
    // Initialize mapped parameters if they were not set
    params.atac.reference   = params.atac.reference ?: [:]

    // Set global assay params
    params.atac.assay                     = 'ATAC'

    GET_BLOCKLIST_FILES(
        params.atac
    )

    validator = new Validate(workflow, params, params.atac, log, messages)
    validator.run()

    // Validate reference workflow
    if (params.atac.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.atac.assay] Executing reference workflow.")
        messages.add("INFO: [$params.atac.assay] Executing reference workflow.")
    }

    // Validate analysis workflow
    if (params.atac.workflow in [ 'analysis', 'full' ]) {
        log.info("INFO: [$params.atac.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.atac.assay] Executing analysis workflow.")
        GET_FASTQ_FILES(
            params.atac
        )
        validator.runAnalysisValidation()
    }

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
        name: "${params.atac.reportsDir}/messages.txt", newLine: true, sort: 'index')

    // Generate references
    if (params.atac.workflow in [ 'reference', 'full' ]) {
        // Set workflow channels
        ch_atac_reference_fasta = Channel.fromPath(params.atac.reference.fasta)
        ch_atac_reference_gtf   = Channel.fromPath(params.atac.reference.gtf)

        // Execute reference workflow
        ATAC_NORMAL_REFERENCE(
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_images_pulled
        )
    }

    // Analyze data
    if (params.atac.workflow in [ 'analysis', 'full' ]) {
        // Set workflow channels
        ch_atac_fastq  = Channel.fromFilePairs("${params.atac.input}" + Core.fastqRegEx(), flat: true)
            .map { prefix, file1, file2 -> tuple(prefix - ~/(_L[0-9][0-9][0-9])$/, file1, file2) }
            .groupTuple()
            .filter { !(it =~ /Undetermined/) }
            .ifEmpty {
                log.error("ERROR: [$params.atac.assay] " + Core.nameFormatMessage())
                exit(1)
            }
        ch_atac_params_file                     = params.atac.paramsFile == null ? Channel.empty()          :
            Channel.from(params.atac.paramsFile)
        if (params.atac.mixed) {
            ch_atac_reference_fasta = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.fasta      :
                Channel.fromPath("${params.atac.reference.directory}/mixed.fa").first()
            ch_atac_reference_gtf   = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.gtf        :
                Channel.fromPath("${params.atac.reference.directory}/mixed.filtered.gtf").first()
        } else {
            if (!(workflow.profile =~ /(awsbatch|tower)/)) { // groovylint-disable-line
                ch_atac_reference_fasta = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.fasta  :
                    Channel.fromPath(params.atac.reference.fasta).first()
                ch_atac_reference_gtf   = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.gtf    :
                    Channel.fromPath(params.atac.reference.gtf).first()
            } else {
                ch_atac_reference_fasta = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.fasta  :
                    Channel.fromPath("${params.atac.reference.directory}/*.{fa,fasta}").first()
                ch_atac_reference_gtf   = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.gtf    :
                    Channel.fromPath("${params.atac.reference.directory}/*.gtf").first()
            }
        }
        if (!(workflow.profile =~ /(awsbatch|tower)/)) { // groovylint-disable-line
            ch_atac_reference_blocklist = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.blocklist  :
            Channel.fromPath(params.atac.reference.blocklist).first()
        } else {
            ch_atac_reference_blocklist = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.blocklist  :
            Channel.fromPath(params.atac.reference.blocklist)
        }
        ch_atac_reference_index     = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.bwaIndex  :
            Channel.fromPath("${params.atac.reference.directory}/bwa-index").first()
        ch_atac_reference_sizes     = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.genomeSize :
            Channel.fromPath("${params.atac.reference.directory}/genome.sizes").first()
        /* groovylint-disable */
        ch_atac_archr_ref           = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.archrRef   :
            Channel.fromPath("${params.atac.reference.directory}/archr/BSgenome.ref.na.1.0_1.0.tar.gz")
        /* groovylint-enable */
        ch_atac_tss                 = params.atac.workflow == 'full' ? ATAC_NORMAL_REFERENCE.out.tss        :
            Channel.fromPath("${params.atac.reference.directory}/TSS.${params.atac.tssWindowSize}.bed")

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
                log.error('ERROR: [ATAC] Missing ATAC reference file(s). Check parameters and/or run ' + // groovylint-disable-line
                    'reference workflow.\n  ' + it.join('\n  ')) // groovylint-disable-line
                exit(1)
            }

        // Execute workflow
        ATAC_NORMAL_ANALYSIS(
            ch_atac_fastq,
            ch_atac_params_file,
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_atac_reference_index,
            ch_atac_reference_sizes,
            ch_atac_archr_ref,
            ch_atac_reference_blocklist,
            ch_atac_tss,
            ch_atac_messages,
            params.atac.assay,
            ch_images_pulled
        )
    }
}
