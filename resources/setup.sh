#!/bin/bash

# Check scripts up to date
git pull

# Initialize variables
config_file="./config"

# Parse options using getopts
while getopts "c:" opt; do
    case $opt in
        c) config_file=$OPTARG ;;
        *) echo "Usage: $0 -c <config_file>"
           exit 1 ;;
    esac
done

# Shift option arguments, so $1 becomes the first positional argument
shift $((OPTIND - 1))

# Initialize an empty string
concatenated=""

# Loop through all arguments
for arg in "$@"; do
    concatenated="$concatenated$arg "
done


set -e
echo "-----------------------------------------------"
echo ""
echo "Using config located at:" ${config_file}
echo ""
echo "-----------------------------------------------"
	
source ${config_file}

if [ -n "${APPTAINER_BIN:-}" ] && [ -n "${HASE_SIF:-}" ]; then
    apptainer_bind="${APPTAINER_BIND:-${home_directory},${scripts_directory}}"
    PYTHON_RUNNER=("${APPTAINER_BIN}" exec)
    if [ -n "${apptainer_bind}" ]; then
        PYTHON_RUNNER+=(--bind "${apptainer_bind}")
    fi
    PYTHON_RUNNER+=("${HASE_SIF}" python)
else
    PYTHON_RUNNER=("${Python_directory}python")
fi
