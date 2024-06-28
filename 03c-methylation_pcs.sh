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
exec &> >(tee ${section_03c_logfile})
print_version

phenfile="NULL"

echo "Calculating methylation PCs of transformed DNA methylation data"
${R_directory}Rscript resources/methylation/methylation_pcs.R \
    ${transformed_methylation_adjusted}.RData \
    ${meth_pc_cutoff} \
    ${n_meth_pcs} \
    ${phenfile} \
    ${meth_pcs_transformed}

echo "Calculating methylation PCs of untransformed DNA methylation data"
${R_directory}Rscript resources/methylation/methylation_pcs.R \
    ${untransformed_methylation_adjusted}.RData \
    ${meth_pc_cutoff} \
    ${n_meth_pcs} \
    ${phenfile} \
    ${meth_pcs_untransformed}

echo "Successfully generated non-genetic methylation PCs"
