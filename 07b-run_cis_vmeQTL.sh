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
exec &> >(tee ${section_07b_logfile}_$1_$2_$3)
print_version

vQTL_method=$1
chr=$2
genetic_chunk=$3

echo "Running cis vmeQTL analysis - CpGs on chr ${chr}; genetic chunk ${genetic_chunk}; method ${vQTL_method}."

if [ $vQTL_method=="BF" ]
then
${osca} \
    --vqtl \
    --vqtl-mtd 2 \
    --bfile ${tabfile}.tab.${genetic_chunk} \
    --befile ${meth_vmeQTL_input_chr}${chr} \
    --cis \
    --cis-wind 2000 \
    --thread-num ${vmeQTL_thread} \
    --task-num 1 \
    --task-id 1 \
    --out ${section_07_dir}/vQTL_${vQTL_method}_cis_genetic${genetic_chunk}_cpgchr${chr}
else
${osca} \
    --vqtl \
    --vqtl-mtd ${vQTL_method} \
    --geno ${tabfile}.tab.${genetic_chunk} \
    --pheno-bod ${meth_vmeQTL_input_chr}${chr} \
    --cis \
    --cis-wind 2000000 \
    --thread-num ${vmeQTL_thread} \
    --task-num 1 \
    --task-id 1 \
    --out ${section_07_dir}/vQTL_${vQTL_method}_cis_genetic${genetic_chunk}_cpgchr${chr}
fi
