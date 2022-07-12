/*
Remove defined number of bases from beginning of R2 reads
*/

params.options = [:]

process CUTADAPT_HEADCROP {
    tag "${sample_id}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    if (workflow.profile == 'aws') {
        label 'large'
    } else {
        label 'cpu_medium'
        label 'memory_medium'
    }

    input:
    tuple val(sample_id), path(fastq)
    val images_pulled

    output:
    tuple val(sample_id), path('*.headcrop.fastq.gz'), emit: fastq
    tuple val(sample_id), path('*_cutadapt_headcrop.log'), emit: log
    tuple val(sample_id), path('*_cutadapt_headcrop_read_counts.csv'), emit: count
    tuple val(sample_id), path('*_cutadapt_headcrop_r2_length.csv'), emit: readlength

    script:
    """
    # Removing bases from 5' end of R2 reads
    cutadapt \
        -U ${params.options.trim.get( sample_id )} \
        -m :20 \
        -o ${sample_id}_R1.headcrop.fastq.gz \
        -p ${sample_id}_R2.headcrop.fastq.gz \
        -j ${task.cpus} \
        ${fastq[0]} \
        ${fastq[1]} \
        1> ${sample_id}_cutadapt_headcrop.log

    # Write average R2 length post-trimming to file
    READLENGTHFILE="${sample_id}_cutadapt_headcrop_r2_length.csv"
    printf "sample,avg_length\n" > \$READLENGTHFILE
    printf "${sample_id},\$(zcat ${sample_id}_R2.headcrop.fastq.gz | awk '{if(NR%4==2) \
        {count++; bases += length} } END{print bases/count}')\n" >> \$READLENGTHFILE

    # Write number of input and output reads to file
    READCOUNTFILE="${sample_id}_cutadapt_headcrop_read_counts.csv"
    printf "sample,process,metric,count\n" > \$READCOUNTFILE
    printf "${sample_id},cutadapt_headcrop,input,\$(grep 'Total read pairs processed:' \
        ${sample_id}_cutadapt_headcrop.log  | awk '{print \$NF}' | sed 's/,//g')\n" \
        >> \$READCOUNTFILE
    printf "${sample_id},cutadapt_headcrop,output,\$(grep 'Pairs written (passing filters):' \
        ${sample_id}_cutadapt_headcrop.log | awk '{print \$5}' | sed 's/,//g')\n" \
        >> \$READCOUNTFILE
    """
}
