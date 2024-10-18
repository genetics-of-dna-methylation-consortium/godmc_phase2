#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_07d_logfile})
print_version

cd $home_directory

suff="tgz"
flags="czf"

tar ${flags} ${home_directory}/results/${study_name}_07.${suff} ${scripts_directory}/config ${scripts_directory}/resources/parameters ${home_directory}/results/07

echo "Successfully created results archives of module 07"
