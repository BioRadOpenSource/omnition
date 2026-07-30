/*
Creating gene symbols files for mixed species runs
*/

process MIXED_GENE_SYMBOLS {
    container "${params.container_image.core_public}:${workflow.manifest.version}"
    publishDir "${params.options.reference.directory}/", mode: 'copy', overwrite: true
    label 'cpu_xsmall'
    label 'memory_xsmall'

    input:
    path gtf
    val images_pulled

    output:
    path '*_gene_symbols.txt', emit: symbols

    script:
    """
    # Creating a text file of the gene symbols for count matrix builds for species 1
    awk 'BEGIN{OFS="\\t";} {if( \$1 ~ /${params.options.species[0]}/) {if ( \$3 ~ /exon|utr/ ) \
        {id=""; name=""; {for(i=1;i<=NF;i++) {if(\$i=="gene_id"){id=\$(i+1)} else if (\$i=="gene_name") \
        {name=\$(i+1)}}} {if(name=="") {print id, id} else {print id, name}}}}}' ${gtf} | tr -d '";' \
        | sort -T '.' | uniq > ${params.options.species[0]}_gene_symbols.txt

    # Creating a text file of the gene symbols for count matrix builds for species 2
    awk 'BEGIN{OFS="\\t";} {if( \$1 ~ /${params.options.species[1]}/) {if ( \$3 ~ /exon|utr/ ) {id=""; \
        name=""; {for(i=1;i<=NF;i++) {if(\$i=="gene_id"){id=\$(i+1)} else if (\$i=="gene_name") {name=\$(i+1)}}} \
        {if(name=="") {print id, id} else {print id, name}}}}}' ${gtf} | tr -d '";' | sort -T '.' | uniq > \
        ${params.options.species[1]}_gene_symbols.txt
    """

    stub:
    """
    # Create expected output files
    outputs=(
        "${params.options.species[0]}_gene_symbols.txt"
        "${params.options.species[1]}_gene_symbols.txt"
    )

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("params.options.species[0]" "params.options.species[1]"
    "gtf" "workflow.manifest.version" "params.options.reference.directory")

    parameters=(
        "${params.options.species[0]}"
        "${params.options.species[1]}"
        "${gtf}"
        "${workflow.manifest.version}"
        "${params.options.reference.directory}"
    )

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "MIXED_GENE_SYMBOLS" "\$param_name" "\$parameter"
    done
    """
}
