#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_02_dir}/logs_c/"

exec &> >(tee ${section_02c_logfile})
print_version

# please install hail from https://hail.is/docs/0.2/getting_started.html before running this script

# Activate the environment
if [ -z "$Python3_directory" ]; then
    if command -v conda &> /dev/null; then
        if conda env list | grep -qE 'hail_env(\s|$)'; then
            echo "found hail_env environment in conda"
            echo "Using conda to run the script"
            RUN_CMD="conda run -n hail_env"
        elif command -v mamba &> /dev/null; then        
            if mamba env list | grep -qE 'hail_env(\s|$)'; then
                echo "found hail_env environment in mamba"
                echo "Using mamba to run the script"
                RUN_CMD="mamba run -n hail_env"
            else
                echo "ERROR: hail_env environment does not exist in conda or mamba. Please ensure hail_env exists."
                exit 1
            fi
        else
            echo "ERROR: hail_env environment does not exist in conda, and mamba not found."
            exit 1
        fi
    elif command -v mamba &> /dev/null; then
        if mamba env list | grep -qE 'hail_env(\s|$)'; then
            echo "found hail_env environment in mamba"
            echo "Using mamba to run the script"
            RUN_CMD="mamba run -n hail_env"
        else
            echo "ERROR: hail_env environment does not exist in mamba. Please ensure hail_env exists."
            exit 1
        fi
    else
        echo "ERROR: Neither conda nor mamba found. Please specify your Python3_directory in the config file."
        exit 1
    fi
else
    echo "Using specified Python3_directory"
    RUN_CMD="${Python3_directory}"
    if ! command -v $RUN_CMD &> /dev/null; then
        echo "ERROR: Specified Python3_directory ($RUN_CMD) is not a valid executable."
        exit 1
    fi
fi

# echo "Running global PCA with raw data"
# echo "RUN_CMD is: $RUN_CMD"
# $RUN_CMD python "${scripts_directory}/resources/datacheck/global_pca.py" \
#     "${section_02_dir}/logs_c/hail_raw.log" \
#     "${bfile_raw}" \
#     "${study_name}" \
#     "${home_directory}" \
#     "${scripts_directory}"

echo "Running global PCA with cleaned data from 02a"
$RUN_CMD python "${scripts_directory}/resources/datacheck/global_pca.py" \
    "${section_02_dir}/logs_c/hail_cleaned.log" \
    "${bfile}" \
    "${study_name}" \
    "${home_directory}" \
    "${scripts_directory}"

echo "Successfully completed godmc2 02c"