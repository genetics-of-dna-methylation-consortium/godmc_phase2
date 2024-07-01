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
exec &> >(tee ${section_06_logfile})
print_version

${R_directory}/Rscript ./resources/methylation/cell_interaction_process.R \
    ${untransformed_methylation_adjusted}.RData \
    ${methylation_array} \
    ${cell_interaction}

for chr in $(seq 1 22)
do
${plink} --bfile ${bfile} \
    --chr $chr \
    --make-bed --out ${tabfile}.chr${chr}
done

