#!/bin/bash

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
mkdir -p ${section_08_dir}/logs_a/
set -e
exec &> >(tee ${section_08a_logfile})
print_version

# Change annotation ####
#'#################################################################################

echo "Adapting plink files"

mkdir -p ${inv_processed_dir}
## Remove indels
${plink2} --bfile  ${bfile} --extract range resources/inversions/inversion_ranges.txt --make-bed --out ${inv_processed_dir}/invSNPs

echo "Inferring inversions"

${R_directory}Rscript --vanilla resources/inversions/genotypeInversions.R ${inv_processed_dir}/invSNPs ${inv_processed_dir} ${section_08_dir} ${nthreads} ${inversion_chunks}

echo "Successfully inferred inversions"


