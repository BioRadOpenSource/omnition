/*
Validate parameters and execute atac combinatorial workflows
*/

// Import workflows
include { PULL_IMAGES } from '../core/singularity.nf'
include { ATAC_COMBINATORIAL_ANALYSIS } from './atac_combinatorial_analysis.nf'
include { ATAC_COMBINATORIAL_REFERENCE } from './atac_combinatorial_reference.nf'
include { GET_FASTQ_FILES } from '../core/fastq_files_array.nf'
include { GET_BLOCKLIST_FILES } from '../core/blocklist_files.nf'

workflow ATAC_COMBINATORIAL {
    take:
    messages

    main:
    // Initialize mapped parameters if they were not set
    params.catac.reference   = params.catac.reference ?: [:]

    // Set global assay params
    params.catac.assay                     = 'cATAC'

    GET_BLOCKLIST_FILES(
        params.catac
    )

    validator = new Validate(workflow, params, params.catac, log, messages)
    validator.run()

    // Validate reference workflow
    if (params.catac.workflow in [ 'reference', 'full' ]) {
        log.info("INFO: [$params.catac.assay] Executing reference workflow.")
        messages.add("INFO: [$params.catac.assay] Executing reference workflow.")
    }

    // Validate analysis workflow
    if (params.catac.workflow in [ 'analysis', 'full' ]) {
        log.info("INFO: [$params.catac.assay] Executing analysis workflow.")
        messages.add("INFO: [$params.catac.assay] Executing analysis workflow.")
        GET_FASTQ_FILES(
            params.catac
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
        name: "${params.catac.reportsDir}/messages.txt", newLine: true, sort: 'index')

    // Generate references
    if (params.catac.workflow in [ 'reference', 'full' ]) {
        // Set workflow channels
        ch_atac_reference_fasta = Channel.fromPath(params.catac.reference.fasta)
        ch_atac_reference_gtf   = Channel.fromPath(params.catac.reference.gtf)

        // Execute reference workflow
        ATAC_COMBINATORIAL_REFERENCE(
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_images_pulled
        )
    }

    // Analyze data
    if (params.catac.workflow in [ 'analysis', 'full' ]) {
        // Set workflow channels
        ch_atac_fastq  = Channel.fromFilePairs("${params.catac.input}" + Core.fastqRegEx(), flat: true)
            .map { prefix, file1, file2 -> tuple(prefix - ~/(_L[0-9][0-9][0-9])$/, file1, file2) }
            .groupTuple()
            .filter { !(it =~ /Undetermined/) }
            .ifEmpty {
                log.error("ERROR: [$params.catac.assay] " + Core.nameFormatMessage())
                exit(1)
            }
        ch_atac_params_file                     = params.catac.paramsFile == null ? Channel.empty()           :
            Channel.from(params.catac.paramsFile)
        ch_atac_TI_config = params.catac.barcodedTn5Config ? Channel.fromPath(params.catac.barcodedTn5Config) :
            Channel.empty()
        if (params.catac.mixed) {
            ch_atac_reference_fasta = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.fasta      :
                Channel.fromPath("${params.catac.reference.directory}/mixed.fa").first()
            ch_atac_reference_gtf   = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.gtf        :
                Channel.fromPath("${params.catac.reference.directory}/mixed.filtered.gtf").first()
        } else {
            if (!(workflow.profile =~ /(awsbatch|tower)/)) { // groovylint-disable-line
                ch_atac_reference_fasta = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.fasta  :
                    Channel.fromPath(params.catac.reference.fasta).first()
                ch_atac_reference_gtf   = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.gtf    :
                    Channel.fromPath(params.catac.reference.gtf).first()
            } else {
                ch_atac_reference_fasta = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.fasta  :
                    Channel.fromPath("${params.catac.reference.directory}/*.{fa,fasta}").first()
                ch_atac_reference_gtf   = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.gtf    :
                    Channel.fromPath("${params.catac.reference.directory}/*.gtf").first()
            }
        }
        if (!(workflow.profile =~ /(awsbatch|tower)/)) { // groovylint-disable-line
            ch_atac_reference_blocklist = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.blocklist :
            Channel.fromPath(params.catac.reference.blocklist).first()
        } else {
            ch_atac_reference_blocklist = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.blocklist :
            Channel.fromPath(params.catac.reference.blocklist)
        }
        ch_atac_reference_index     = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.bwaIndex  :
            Channel.fromPath("${params.catac.reference.directory}/bwa-index").first()
        ch_atac_reference_sizes     = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.genomeSize :
            Channel.fromPath("${params.catac.reference.directory}/genome.sizes").first()
        /* groovylint-disable */
        ch_atac_archr_ref           = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.archrRef   :
            Channel.fromPath("${params.catac.reference.directory}/archr/BSgenome.ref.na.1.0_1.0.tar.gz")
        /* groovylint-enable */
        ch_atac_tss                 = params.catac.workflow == 'full' ? ATAC_COMBINATORIAL_REFERENCE.out.tss        :
            Channel.fromPath("${params.catac.reference.directory}/TSS.${params.catac.tssWindowSize}.bed")

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

        ATAC_COMBINATORIAL_ANALYSIS(
            ch_atac_fastq,
            ch_atac_params_file,
            ch_atac_reference_fasta,
            ch_atac_reference_gtf,
            ch_atac_reference_index,
            ch_atac_reference_sizes,
            ch_atac_archr_ref,
            ch_atac_reference_blocklist,
            ch_atac_tss,
            ch_atac_TI_config,
            ch_atac_messages,
            params.catac.assay,
            ch_images_pulled
        )
    }
}
