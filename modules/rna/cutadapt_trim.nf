/*
Trimming adapters, low-quality regions, and polyA sequences
*/

params.options = [:]

process CUTADAPT_TRIM {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(fastq), path(fastq_valid)
    val images_pulled

    output:
    tuple val(sampleId), path('*.trimmed.fastq.gz'), emit: fastq
    tuple val(sampleId), path('*_cutadapt_trim_read_counts.csv'), emit: count
    tuple val(sampleId), path('*_cutadapt_trim.log'), emit: log
    tuple val(sampleId), path('*_cutadapt_trim_r2_length.csv'), optional: true, emit: readlength

    script:
    """
    # Trimming reads, none of these seqs should be here but have nonetheless been observed
    cutadapt \
        -G "Trueseq_P5_Forward=ACGACGCTCTTCCGATCT;max_errors=0;min_overlap=5" \
        -m :20 \
        -q ${params.readQualityScore} \
        -o ${sampleId}_R1.trimmed.fastq.gz \
        -p ${sampleId}_R2.trimmed.fastq.gz \
        -n 4 \
        -j ${task.cpus} \
        ${fastq} \
        1> ${sampleId}_cutadapt_trim.log

    # Write average R2 length post-trimming to file
    READLENGTHFILE="${sampleId}_cutadapt_trim_r2_length.csv"
    printf "sample,avg_length\n" > \$READLENGTHFILE
    printf "${sampleId},\$(zcat ${sampleId}_R2.trimmed.fastq.gz | awk '{if(NR%4==2) \
        {count++; bases += length} } END{print bases/count}')\n" >> \$READLENGTHFILE

    # Write number of input and output reads to file
    READCOUNTFILE="${sampleId}_cutadapt_trim_read_counts.csv"
    printf "sample,process,metric,value\n" > \$READCOUNTFILE
    printf "${sampleId},cutadapt_trim,input,\
        \$(grep 'Total read pairs processed:' ${sampleId}_cutadapt_trim.log \
        | awk '{print \$NF}' | sed 's/,//g')\n" >> \$READCOUNTFILE
    printf "${sampleId},cutadapt_trim,output,\
        \$(grep 'Pairs written (passing filters):' ${sampleId}_cutadapt_trim.log \
        | awk '{print \$5}' | sed 's/,//g')\n" >> \$READCOUNTFILE
    """
}
