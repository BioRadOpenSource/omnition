/*
Filtering GTF entries to only contain chromosomes/contigs in the FASTA file
*/

params.options = [:]

process FILTER_REFERENCES {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_small'

    input:
    path fasta
    path gtf
    val images_pulled

    output:
    path "${gtf.baseName}.filtered.gtf", emit: gtf

    script:
    """
    # Transfering over all commented lines into the new filtered file
    awk '/^#/' ${gtf} > ${gtf.baseName}.filtered.gtf

    # Check if the input files are empty
    if  ! [ -s ${fasta} ] && [ -s ${gtf} ]; then
        echo "[ERROR] gtf or fasta input file empty before filtering"
        exit 1
    else
        echo "data found in both fasta and gtf files"
    fi

    # Creates an array from lines in the first file that start with > and puts into it
    # all the values in the first element after the >.
    # Then for only the second file it looks to see if the first element in each line of the
    # second file is in the array previously created and only transfers those into the filtered file
    awk '/>/ && NR==FNR{ array[substr(\$1, 2, length(\$1))]; next } \$1 in array' \
        ${fasta} ${gtf} >> ${gtf.baseName}.filtered.gtf

    # Check if the filtered gtf is completely empty
    if ! [ -s ${gtf.baseName}.filtered.gtf ]; then
        echo "[ERROR] filtered gtf file empty"
        exit 1
    fi
    """
}
