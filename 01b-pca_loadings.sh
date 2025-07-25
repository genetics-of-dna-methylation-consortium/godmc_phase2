#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_01_dir}/logs_b/"

exec &> >(tee ${section_01b_logfile})
print_version

# please install hail from https://hail.is/docs/0.2/getting_started.html before running this script

# Activate the environment
if [ -z "$Python3_directory" ]; then
    # Users are using system default R & Python, activate conda environment
    if command -v mamba &> /dev/null; then
        CONDA_CMD="mamba"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
    else
        echo "Neither mamba nor conda found. Please install one of them first."
        exit 1
    fi
    # Initialize shell for current session
    eval "$($CONDA_CMD shell hook --shell bash)"
    $CONDA_CMD activate hail_env
    echo "Current conda environment: $CONDA_DEFAULT_ENV"
else
    # Users have specified custom R/Python directories
    echo "Custom R/Python directories specified"
fi

echo "Running global PCA"
${Python3_directory}python "${scripts_directory}/resources/datacheck/global_pca.py" \
    "${section_01_dir}/logs_b/hail.log" \
    "${bfile_raw}" \
    "${study_name}" \
    "${scripts_directory}"

echo "Successfully completed godmc2 01b"