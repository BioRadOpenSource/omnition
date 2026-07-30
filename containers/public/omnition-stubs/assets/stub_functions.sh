#!/bin/bash

# These are a set of functions that are used by stub blocks to:
# Record parameters used in module (check_params)
# Check if expected # of input files match with what's expected (check_input_files)

function check_params() {
    MODULE=$1
    PARAMETER=$2
    VALUE=$3
    # Collect all parameters and values
    echo "$PARAMETER: $VALUE" >> parameter_stubs.txt

    # Give warning if input file/parameter is missing
    if [[ ! -f "$VALUE" && -z "$VALUE" ]]; then
        echo "The following input file/parameter is missing:";
        echo "Param/file: $PARAMETER with value: $VALUE";
        # Collect warnings with module, variable, and value
        echo "$MODULE $PARAMETER $VALUE" >> warning_stub.txt
    fi;
}; export -f check_params

function check_input_files() {
    MODULE=$1
    sampleId=$2
    EXPECTED_INPUT_COUNT=$3

    count=$(find . -type f -name "*${sampleId}*" | wc -l)

    if [ "$count" -ne "$EXPECTED_INPUT_COUNT" ]; then
        echo "$MODULE Error: $sampleId has $count input files,
        It is specified to have $EXPECTED_INPUT_COUNT input files." \
        >> error_stub.txt
    fi
}; export -f check_input_files
