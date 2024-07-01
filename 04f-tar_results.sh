#!/bin/bash -l

set -e

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

set -e
echo "-----------------------------------------------"
echo ""
echo "Using config located at:" ${config_file}
echo ""
echo "-----------------------------------------------"

source ${config_file}

exec &> >(tee ${section_04f_logfile})
print_version

cd $home_directory

suff="tgz"
flags="czf"

tar ${flags} ${home_directory}/results/${study_name}_04.${suff} ${scripts_directory}/config ${scripts_directory}/resources/parameters ${home_directory}/results/04

echo "Successfully created results archives of module 04"
