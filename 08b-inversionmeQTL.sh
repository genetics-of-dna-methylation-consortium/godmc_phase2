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
mkdir -p ${section_08_dir}
mkdir -p ${section_08_dir}/logs_b/
set -e
exec &> >(tee ${section_08b_logfile})
print_version

geno="${inv_processed_dir}/inversionQTL.txt"
methy="${transformed_methylation_adjusted_pcs}.RData"
methy_trans="${inv_processed_dir}/methy_transposed.txt"

out="${section_08_dir}/invmeqtl.Rdata"
threshold=1

echo "Transposing methylation matrix for MatrixEQTL"
Rscript resources/inversions/adapt_methylation_meQTL.R ${methy} ${methy_trans}

echo "Performing inversion meQTL analysis"

${R_directory}Rscript resources/genetics/matrixeqtl.R ${geno} ${methy_trans} NULL ${threshold} ${out}


echo "Successfully completed inversion meQTL analysis"
