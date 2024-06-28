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
exec &> >(tee ${section_03f_logfile})
print_version

mkdir -p ${hase_cov}
mkdir -p ${hase_cov_males}
mkdir -p ${hase_cov_females}

${R_directory}Rscript resources/methylation/methylation_hase_format.R \
		${untransformed_methylation_adjusted_pcs} \
		${transformed_methylation_adjusted_pcs} \
		${covariates_combined}.txt \
		${hase_cov} \
                ${hase_cov_males} \
                ${hase_cov_females}

echo "Successfully converted methylation data"

