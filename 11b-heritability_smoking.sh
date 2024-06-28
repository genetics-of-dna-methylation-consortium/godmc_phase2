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
mkdir -p ${section_11_dir}/logs_b
touch ${section_11b_logfile}
exec &> >(tee ${section_11b_logfile})
print_version

#### GREML for SNP heritability###################################

  ${gcta} \
          --grm ${grmfile_all}  \
          --reml \
          --pheno ${smoking_pred}.smok.plink \
          --out ${section_11_dir}/heritability_smoking \
          --thread-num ${nthreads}
          
echo "Successfully finished the calculation on SNP heritability for smoking!"          
        