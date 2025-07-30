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
        # Initialize shell for current session
        eval "$($CONDA_CMD shell hook --shell bash)"
    elif command -v conda &> /dev/null; then
        CONDA_CMD="conda"
        # Initialize shell for current session
        eval "$($CONDA_CMD shell.bash hook)"
    else
        echo "ERROR: Neither mamba nor conda found."
        echo "Please install one of them first or specify the R/Python directories in the config file"
        exit 1
    fi

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
    "${home_directory}"

echo "Successfully completed godmc2 01b"