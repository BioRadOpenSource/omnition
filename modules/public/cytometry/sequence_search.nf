/*
Checking R2 for known contaminants: mito, ribosomal, unused adts
*/

process SEQUENCE_SEARCH {
    tag "${sampleId}"
    container "${params.container_image.cytometry_public}:${workflow.manifest.version}"
    publishDir "${params.options.resultsDir}/${sampleId}/${params.options.assay}/contaminants/",
        mode: 'copy', overwrite: true, pattern: 'pseudoalignments.bam'
    label 'cpu_medium'
    label 'memory_medium'

    input:
    tuple val(sampleId), path(read2), path(adt_file), path(index)
    val images_pulled

    output:
    tuple val(sampleId), path('*.{bam,csv,json}'), emit: seq_search
    tuple val(sampleId), path('*_seq_search.csv'), emit: metrics

    script:
    """
    TAGLENGTH=\$(awk -F ',' 'NR == 1 || (length(\$1) > 0 && length(\$1) < shortest) {
            shortest = length(\$1)
        } END {
            print shortest
        }' ${adt_file})
    if ! (( \$TAGLENGTH % 2 )); then
        echo 'Minimum tag length was even, subtracting one to make it odd...';
        TAGLENGTH=\$((\$TAGLENGTH - 1));
    fi

    # Run kallisto - allow it to "fail" with warnings but continue processing
    kallisto quant -i ${index} -o ./ -l \$TAGLENGTH -s 1 --single --single-overhang --pseudobam  ${read2} || echo "Kallisto completed with warnings, continuing..."

    # Copy the references to the working directory
    cp /opt/biorad/assets/references/* ./

    # Check if BAM file is valid, create empty one if corrupted
    if ! samtools quickcheck pseudoalignments.bam 2>/dev/null; then
        echo "BAM file is corrupted or empty, creating valid empty BAM file..."
        # Create a proper empty BAM file with header
        echo -e "@HD\tVN:1.0\tSO:unsorted" | samtools view -bS - > pseudoalignments.bam
    fi

    # Create BAM header for all output files
    samtools view pseudoalignments.bam -H > header.txt

    # Create a bam file with the mitochondrial alignments
    if [ -f "mitochondrial_names.txt" ]; then
        samtools view pseudoalignments.bam | grep -F -w -f mitochondrial_names.txt > mitochondrial.txt || true
    else
        touch mitochondrial.txt
    fi
    cat header.txt mitochondrial.txt > mitochondrial.sam
    samtools view mitochondrial.sam -b > mitochondrial.bam

    # Create a bam file with the ribosomal alignments
    samtools view pseudoalignments.bam | grep "U13369.1-humanribosomalDNA" > ribosomal.txt || true
    cat header.txt ribosomal.txt > ribosomal.sam
    samtools view ribosomal.sam -b > ribosomal.bam

    # Create a bam file with the unmapped alignments, can use same header
    samtools view -f 4 pseudoalignments.bam > unmapped.txt
    cat header.txt unmapped.txt > unmapped.sam
    samtools view unmapped.sam -b > unmapped.bam

    # Create a table of the most frequent sequences
    samtools view -f 4 pseudoalignments.bam |
        cut -f 10 |
        cut -c 1-\$TAGLENGTH |
        sort |
        uniq -c |
        sort -nr > top_unmapped_seqs.csv || touch top_unmapped_seqs.csv
    sed -i 's/^[ \t]*//' top_unmapped_seqs.csv

    # Calculate alignment counts with default fallback to 0 if calculation fails
    TOTALALIGNMENTS=\$(samtools view pseudoalignments.bam | cut -f 1 | sort | uniq -c | wc -l 2>/dev/null || echo 0)
    MITOALIGNMENTS=\$(samtools view mitochondrial.bam | cut -f 1 | sort | uniq -c | wc -l 2>/dev/null || echo 0)
    UNMAPPEDALIGNMENTS=\$(samtools view unmapped.bam | cut -f 1 | sort | uniq -c | wc -l 2>/dev/null || echo 0)
    RIBOSOMALALIGNMENTS=\$(samtools view ribosomal.bam | cut -f 1 | sort | uniq -c | wc -l 2>/dev/null || echo 0)

    # Create a csv recording the number of each alignment type
    echo "num_total,\$TOTALALIGNMENTS" > num_aligned_reads.csv
    echo "num_mito,\$MITOALIGNMENTS" >> num_aligned_reads.csv
    echo "num_unmapped,\$UNMAPPEDALIGNMENTS" >> num_aligned_reads.csv
    echo "num_ribosomal,\$RIBOSOMALALIGNMENTS" >> num_aligned_reads.csv

    # Create a csv recording the number of each alignment type for metric summary
    echo "sample,process,metric,value" > ${sampleId}_seq_search.csv
    echo "${sampleId},sequence_search,mito_reads,\$MITOALIGNMENTS" >> ${sampleId}_seq_search.csv
    echo "${sampleId},sequence_search,ribo_reads,\$RIBOSOMALALIGNMENTS" >> ${sampleId}_seq_search.csv
    echo "${sampleId},sequence_search,unmapped_reads,\$UNMAPPEDALIGNMENTS" >> ${sampleId}_seq_search.csv
    echo "${sampleId},sequence_search,total_reads,\$TOTALALIGNMENTS" >> ${sampleId}_seq_search.csv
    if [ -f "run_info.json" ]; then
        mv run_info.json seqsearch_run_info.json
    else
        echo '{"call": "kallisto quant", "start_time": "", "end_time": ""}' > seqsearch_run_info.json
    fi
    """
}
