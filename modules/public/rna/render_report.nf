/*
Generate report
*/

process RENDER_REPORT {
    container "${params.container_image.customer_report_public}:${workflow.manifest.version}"
    publishDir "${params.options.reportsDir}", mode: 'copy', overwrite: true
    label 'cpu_small'
    label 'memory_medium'

    input:
    path samples
    path settings
    path input
    val images_pulled

    output:
    path("${params.options.prefix}rna_*.html")

    script:
    """
    # If a file starting with metric_summary_ exists, copy its contents to metric_summary.csv
    if ls metric_summary_* 1> /dev/null 2>&1; then
        cat metric_summary_* > metric_summary.csv
    fi

    mkdir tmp
    mkdir -p results/sampleId/plots
    echo "${params}" | sed 's/,/\\n/g' | sed 's/:/,/g' | tr -d '[' | tr -d ']' | tail -n +3 > ./tmp/params.csv

    # If there are any other message files, combine them
    if [ -n "\$(ls *_messages.txt 2>/dev/null)" ];
    then
        cat messages.txt *_messages.txt |
        grep -v '^\$' > "messages.temp" && mv "messages.temp" "messages.txt"
    fi

    # If using Singularity, copy report files to the working directory
    if [[ "${workflow.profile}" =~ "standard" ]]
    then
        # Environment variables are now set in Dockerfile, but adding here as backup
        export MPLCONFIGDIR=\${MPLCONFIGDIR:-/tmp/matplotlib}
        export NPM_CONFIG_CACHE=\${NPM_CONFIG_CACHE:-/tmp/npm_cache}
        export NPM_CONFIG_TMP=\${NPM_CONFIG_TMP:-/tmp/npm_tmp}
        cp -r /opt/biorad/app/* ./
        APP_PREFIX="."
    else
        APP_PREFIX="/opt/biorad/app"
    fi

    # Creating the file that contains the information necessary to display release version
    printf %s "${workflow.manifest.version}" >> version.txt

    # Perform preprocessing of data
    python \${APP_PREFIX}/processing/process_data.py -i . -o \${APP_PREFIX}/app/assets/data/ -t rna

    # If using Docker or Tower, cd to the app directory
    if [[ ${workflow.profile} =~ (docker|tower|awsbatch) ]]
    then
        BASE_DIR=\$PWD
        cd \${APP_PREFIX}
    fi

    # Build the report
    npm run build:rna -- --output-directory=. --filename=${params.options.prefix}rna_\$(date "+%y%m%d-%H%M").html

    # If using Docker or Tower, move the report to the working directory
    if [[ ${workflow.profile} =~ (docker|tower|awsbatch) ]]
    then
        mv \${APP_PREFIX}/*.html \${BASE_DIR}
    fi
    """

    stub:
    """
    # If a file starting with metric_summary_ exists, copy its contents to metric_summary.csv
    if ls metric_summary_* 1> /dev/null 2>&1; then
        cat metric_summary_* > metric_summary.csv
    fi

    # Create expected output files
    outputs=("${params.options.prefix}rna_*.html")

    for output in "\${outputs[@]}"
    do
        touch \$output
    done

    # Record all groovy parameters used in module
    param_names=("workflow.manifest.version" "params.options.prefix"
    "params.options.reportsDir" "workflow.profile" "params")

    parameters=("${workflow.manifest.version}" "${params.options.prefix}"
    "${params.options.reportsDir}" "${workflow.profile}" "${params}")

    # Check that expected params and files exists
    for index in "\${!parameters[@]}"
    do
        param_name="\${param_names[index]}"
        parameter="\${parameters[index]}"
        check_params "RENDER_CUSTOMER_REPORT" "\$param_name" "\$parameter"
    done
    """
}
