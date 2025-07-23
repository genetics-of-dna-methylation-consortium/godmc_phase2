#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_01_dir}/logs_b/"

exec &> >(tee ${section_01b_logfile})
print_version

# please install hail from https://hail.is/docs/0.2/getting_started.html before running this script

mamba activate hail_env

echo "Running global PCA"
${Python3_directory}python "${scripts_directory}/resources/datacheck/global_pca.py" \
    "${section_01_dir}/logs_b/hail.log" \
    "${bfile_raw}" \
    "${study_name}" \
    "${scripts_directory}"

# local pca plot to be done
echo "Successfully completed godmc2 01b"