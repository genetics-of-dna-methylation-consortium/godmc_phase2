#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_16_dir}/logs_g"
exec &> >(tee "${section_16g_logfile}")
print_version

cd "${home_directory}"

suff="tgz"
flags="czf"
archive="${home_directory}/results/${study_name}_16.${suff}"

tar ${flags} "${archive}" \
    --exclude="results/16/positive_control_validation/hase/*/run/*" \
    "${scripts_directory}/config" \
    "${scripts_directory}/resources/parameters" \
    "results/16"

echo "Successfully created results archives of module 16"
