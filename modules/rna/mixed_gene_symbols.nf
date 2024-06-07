/*
Creating gene symbols files for mixed species runs
*/

params.options = [:]

process MIXED_GENE_SYMBOLS {
    container "bioraddbg/omnition-core:${workflow.manifest.version}"
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
}
