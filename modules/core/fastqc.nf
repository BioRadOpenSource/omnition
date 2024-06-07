/*
Assessing read quality
*/

params.options = [:]

process FASTQC {
    tag "${sampleId}"
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    label 'cpu_large'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(fastq), path(fastq_valid)
    val images_pulled

    output:
    tuple val(sampleId), path('*_fastqc.zip'), emit: zip
    tuple val(sampleId), path('*_read_counts.csv'), optional: true, emit: counts
    tuple val(sampleId), path('*_quality_scores.csv'), optional: true, emit: quality_scores
    tuple val(sampleId), path('*_sequence_traces.csv'), optional: true, emit: sequence_traces
    tuple val(sampleId), path('*.html'), optional: true, emit: qc_html

    script:
    // Define common variables
    """
    # Reducing heap size to 80% of allocated resources account for Java overhead
    JAVA_MEMORY="\$(((${task.memory.toGiga()} * 4)/ 5))g"

    if [ ${params.options.fastqcDs} != 1 ]; then
        # Sample FASTQ files to desired read fraction
        reformat.sh \
            in=${fastq[0]} \
            in2=${fastq[1]} \
            out=${sampleId}_R1_001.ds.complete.fastq.gz \
            out2=${sampleId}_R2_001.ds.complete.fastq.gz \
            sampleseed=1 \
            t=${task.cpus} \
            -Xmx\$JAVA_MEMORY \
            2>&1 | tee reformat.log

        # Generate read stats
        fastqc *.ds.complete.fastq.gz --threads ${task.cpus} --quiet  2>&1

        # Extract input read counts to shell variable
        READCOUNTS=`awk '/Input:/{print \$2/2}' reformat.log`
    else
        # Generate read stats
        fastqc ${fastq} --threads ${task.cpus} --quiet  2>&1

        # Extract input read counts to shell variable
        READCOUNTS=\$(grep -oP \
        '(?<=Total Sequences</td><td>).*(?=</td></tr><tr><td>Sequences flag)' *R1_001*_fastqc.html)
    fi

    function process_data {
        # This function pulls together metrics from fastqc into single file
        sampleId=\$1
        readNum=\$2
        COUNTFILE=\$3

        # Convert readNum to lowercase
        readNumLower=\$(echo \$readNum | tr '[:upper:]' '[:lower:]')

        # Extract quality scores to file, taking the last number of a range
        unzip -p "\${sampleId}*_\${readNum}_001*_fastqc.zip" "${sampleId}*_\${readNum}_001*_fastqc/fastqc_data.txt" \
            > \${readNum}_data.txt
        cat \${readNum}_data.txt | awk '/Per base sequence quality/, />>END_MODULE/' | head -n -1 | tail -n +2 | \
            cut -f 1,2 | sed 's;.*-;;' > \${sampleId}_\${readNum}_quality_scores.csv

        # Extract sequence traces to file
        cat \${readNum}_data.txt | awk '/Per base sequence content/, />>END_MODULE/' | head -n -1 | tail -n +2 | \
            sed 's;.*-;;' > \${sampleId}_\${readNum}_sequence_traces.csv

        # Calculate average q30
        ave_q30=\$(cat \${readNum}_data.txt | awk '/Per base sequence quality/, />>END_MODULE/' | head -n -1 | \
            tail -n +3 | cut -f 2 | awk 'NR == 1 { sum=0} {sum+=\$1;} END {printf sum/NR}' )

        # Update read counts file with average q30
        printf "\${sampleId},fastqc,\${readNumLower},q30_\${readNumLower},\$(echo \$ave_q30)\n" >> "\${COUNTFILE}"

        # Write input read counts to file
        printf "${sampleId},fastqc,\${readNumLower},read_count,\$(echo \$READCOUNTS)\n" >> "\${COUNTFILE}"
    }

    # Common code for both if and else blocks
    COUNTFILE=${sampleId}_input_read_counts.csv
    printf "sample,process,read,metric,value\n" > \$COUNTFILE

    # Gather read metrics to file
    process_data "${sampleId}" "R1" "\${COUNTFILE}"
    process_data "${sampleId}" "R2" "\${COUNTFILE}"
    """
}
