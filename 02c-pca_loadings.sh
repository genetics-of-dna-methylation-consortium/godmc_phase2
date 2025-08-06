#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_02_dir}/logs_c/"

exec &> >(tee ${section_02c_logfile})
print_version

# please install hail from https://hail.is/docs/0.2/getting_started.html before running this script

# Activate the environment
if [ -z "$Python3_directory" ]; then
    if command -v mamba &> /dev/null; then
        echo "Using mamba to run the script"
        RUN_CMD="mamba run -n hail_env"
    elif command -v conda &> /dev/null; then
        echo "Using conda to run the script"
        RUN_CMD="conda run -n hail_env"
    else
        echo "ERROR: Neither mamba nor conda found."
        exit 1
    fi
else
    RUN_CMD="${Python3_directory}"
fi

echo "Running global PCA with raw data"
echo "RUN_CMD is: $RUN_CMD"
$RUN_CMD python "${scripts_directory}/resources/datacheck/global_pca.py" \
    "${section_02_dir}/logs_c/hail_raw.log" \
    "${bfile_raw}" \
    "${study_name}" \
    "${home_directory}" \
    "${scripts_directory}"

echo "Running global PCA with cleaned data"
$RUN_CMD python "${scripts_directory}/resources/datacheck/global_pca.py" \
    "${section_02_dir}/logs_c/hail_cleaned.log" \
    "${bfile}" \
    "${study_name}" \
    "${home_directory}" \
    "${scripts_directory}"

echo "Successfully completed godmc2 02c"