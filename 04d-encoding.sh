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
exec &> >(tee ${section_04d_logfile})
print_version

#Please read resources/bin/hase/README_2.md
#An example is also provided below

mkdir -p ${hase_encoding}
mkdir -p ${hase_pheno}
cp ${transformed_methylation_adjusted_pcs}.csv ${hase_pheno}
mv ${hase_pheno}/transformed_methylation_adjusted_pcs.csv ${hase_pheno}/methylation_data.csv

python ${hase}/hase.py \
   -mode encoding \
   -study_name ${study_name} \
   -g ${hase_converting} \
   -o ${hase_encoding} \
   -mapper ${hase_mapping} \
   -ph ${hase_pheno} \
   -ref_name ref-hrc

echo "Successfully encoded the genetic data"
