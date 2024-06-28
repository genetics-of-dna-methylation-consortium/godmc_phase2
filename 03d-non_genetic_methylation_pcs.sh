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
exec &> >(tee ${section_03d_logfile})

print_version

echo "Performing genetic analysis of methylation PCs of transformed data by MatrixEQTL"
${R_directory}Rscript resources/methylation/genetic_meth_pcs.R \
	${tabfile}.tab \
	${meth_pcs_transformed} \
	${nongenetic_meth_pcs_transformed} \
	${meth_pc_threshold} \
	${nthreads}
    
echo "Performing genetic analysis of methylation PCs of untransformed data by MatrixEQTL"
${R_directory}Rscript resources/methylation/genetic_meth_pcs.R \
	${tabfile}.tab \
	${meth_pcs_untransformed} \
	${nongenetic_meth_pcs_untransformed} \
	${meth_pc_threshold} \
	${nthreads}

echo "Successfully generating non-genetic methylation PCs"
