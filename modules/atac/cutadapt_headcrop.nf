/*
Remove defined number of bases from beginning of R2 reads
*/

params.options = [:]

process CUTADAPT_HEADCROP {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(fastq)
    val images_pulled

    output:
    tuple val(sampleId), path('*.headcrop.fastq.gz'), emit: fastq
    tuple val(sampleId), path('*_cutadapt_headcrop.log'), emit: log
    tuple val(sampleId), path('*_cutadapt_headcrop_read_counts.csv'), emit: count
    tuple val(sampleId), path('*_cutadapt_headcrop_r2_length.csv'), emit: readlength

    script:
    """
    # Removing bases from 5' and 3' end of R2 reads
    cutadapt \
        -U ${params.options.trimFivePrime.get( sampleId )} \
        -U -${params.options.trimThreePrime.get( sampleId )} \
        -m :20 \
        -o ${sampleId}_R1.headcrop.fastq.gz \
        -p ${sampleId}_R2.headcrop.fastq.gz \
        -j ${task.cpus} \
        ${fastq[0]} \
        ${fastq[1]} \
        1> ${sampleId}_cutadapt_headcrop.log

    # Write average R2 length post-trimming to file
    READLENGTHFILE="${sampleId}_cutadapt_headcrop_r2_length.csv"
    printf "sample,avg_length\n" > \$READLENGTHFILE
    printf "${sampleId},\$(zcat ${sampleId}_R2.headcrop.fastq.gz | awk '{if(NR%4==2) \
        {count++; bases += length} } END{print bases/count}')\n" >> \$READLENGTHFILE

    # Write number of input and output reads to file
    READCOUNTFILE="${sampleId}_cutadapt_headcrop_read_counts.csv"
    printf "sample,process,metric,value\n" > \$READCOUNTFILE
    printf "${sampleId},cutadapt_headcrop,input,\$(grep 'Total read pairs processed:' \
        ${sampleId}_cutadapt_headcrop.log  | awk '{print \$NF}' | sed 's/,//g')\n" \
        >> \$READCOUNTFILE
    printf "${sampleId},cutadapt_headcrop,output,\$(grep 'Pairs written (passing filters):' \
        ${sampleId}_cutadapt_headcrop.log | awk '{print \$5}' | sed 's/,//g')\n" \
        >> \$READCOUNTFILE
    """
}
