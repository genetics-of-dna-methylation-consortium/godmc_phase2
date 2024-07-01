#!/bin/bash

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
chr=${1}
cell_type=${2}
exec &> >(tee ${section_06b_logfile}_${batch_number})
print_version

${Python_directory}/python resources/methylation/cell_interaction.py \
    ${tabfile} \
    ${cell_interaction} \
    ${cellcounts_cov} \
    ${chr} \
    ${cell_type} \
    ${section_06_dir}
