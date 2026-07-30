/*
Trimming adapters, low-quality regions, and polyA sequences
*/

process CUTADAPT_TRIM {
    tag "${sampleId}"
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    label 'cpu_xsmall'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(fastq), path(fastq_valid)
    val images_pulled

    output:
    tuple val(sampleId), path('*.trimmed.fastq.gz'), emit: fastq
    tuple val(sampleId), path('*_cutadapt_trim_read_counts.csv'), emit: count
    tuple val(sampleId), path('*_cutadapt_trim.log'), emit: log
    tuple val(sampleId), path('*_cutadapt_trim_read_length.csv'), optional: true, emit: readlength

    script:
    """
    # Trimming reads, none of these seqs should be here but have nonetheless been observed
    cutadapt \
        -G "TruSeq_P5_Forward=ACGACGCTCTTCCGATCT;max_errors=0;min_overlap=5" \
        -m :20 \
        -q ${params.readQualityScore} \
        -o ${sampleId}_R1.trimmed.fastq.gz \
        -p ${sampleId}_R2.trimmed.fastq.gz \
        -n 4 \
        -j ${task.cpus} \
        -Z \
        ${fastq} \
        1> ${sampleId}_cutadapt_trim.log

    # Write average read length post-trimming to file
    READLENGTHFILE="${sampleId}_cutadapt_trim_read_length.csv"
    printf "sample,process,metric,value\n" > \$READLENGTHFILE
    printf "${sampleId},cutadapt_trim,avg_length_r1,\$(zcat ${sampleId}_R1.trimmed.fastq.gz | awk '{if(NR%4==2) \
        {count++; bases += length} } END{print bases/count}')\n" >> \$READLENGTHFILE
    printf "${sampleId},cutadapt_trim,avg_length_r2,\$(zcat ${sampleId}_R2.trimmed.fastq.gz | awk '{if(NR%4==2) \
        {count++; bases += length} } END{print bases/count}')\n" >> \$READLENGTHFILE

    # Parse the log file for trimming counts
    POLARS_MAX_THREADS=${task.cpus}
    publicCoreParseCutadaptLog.py -i ${sampleId}_cutadapt_trim.log -s ${sampleId}
    """

    stub:
    """
    # Expected number of sample input files
    EXPECTED_INPUT_COUNT=8

    # Check that we have expected number of sample input files
    check_input_files "CUTADAPT_TRIM" "${sampleId}" "\$EXPECTED_INPUT_COUNT"

    # Create expected output files
    outputs=("${sampleId}_R1.trimmed.fastq.gz" "${sampleId}_R2.trimmed.fastq.gz"
    "${sampleId}_cutadapt_trim_read_counts.csv" "${sampleId}_cutadapt_trim.log"
    "${sampleId}_cutadapt_trim_read_length.csv")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("task.cpus" "params.readQualityScore"
    "workflow.manifest.version"
    "fastq" "params.options.resultsDir" "sampleId")

    parameters=("${task.cpus}"
    "${params.readQualityScore}" "${workflow.manifest.version}"
    "${fastq}"
    "${params.options.resultsDir}" "${sampleId}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "CUTADAPT_TRIM" "\$param_name" "\$parameter"
    done
    """
}
