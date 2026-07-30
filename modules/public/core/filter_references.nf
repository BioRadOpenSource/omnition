/*
Filtering GTF entries to only contain chromosomes/contigs in the FASTA file
*/

process FILTER_REFERENCES {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
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

    stub:
    """
    # Create expected output files
    outputs=("${gtf.baseName}.filtered.gtf")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.reference.directory" "gtf.baseName" "fasta"
    "workflow.manifest.version" "gtf")

    parameters=("${params.options.reference.directory}" "${gtf.baseName}" "${fasta}"
    "${workflow.manifest.version}" "${gtf}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "FILTER_REFERENCES" "\$param_name" "\$parameter"
    done
    """
}
