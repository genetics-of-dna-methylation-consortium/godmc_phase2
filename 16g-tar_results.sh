#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_16_dir}/logs_g"
exec &> >(tee ${section_16g_logfile})
print_version

cd $home_directory

suff="tgz"
flags="czf"

tar ${flags} ${home_directory}/results/${study_name}_16.${suff} ${scripts_directory}/config ${scripts_directory}/resources/parameters ${home_directory}/results/16

echo "Successfully created results archives of module 16"
