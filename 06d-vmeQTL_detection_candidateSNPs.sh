#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}/logs_d
mkdir -p ${section_06_dir}/vmeQTL_results/Trans_candidateSNPs/

exec &> >(tee ${section_06d_logfile}_method${1}_chr${2})
print_version

vQTL_method=$1
chr=$2

echo "Started to run 06d chunk${chunk} method ${vQTL_method} at $(date)"
if [ $vQTL_method = "BF" ]
then
${osca} \
    --vqtl \
    --vqtl-mtd 2 \
    --bfile ${tabfile}.vQTLs \
    --befile ${meth_vmeQTL_input_chr}${chr} \
    --trans \
    --trans-wind 2000 \
    --thread-num 10 \
    --task-num 1 \
    --task-id 1 \
    --out ${section_06_dir}/vmeQTL_results/Trans_candidateSNPs/vQTL_${vQTL_method}_trans_methylation_chr${chr}
else
${osca_new} \
    --vqtl \
    --vqtl-mtd ${vQTL_method} \
    --geno ${tabfile}.vQTLs \
    --pheno-bod ${meth_vmeQTL_input_chr}${chr} \
    --trans \
    --trans-wind 2000000 \
    --thread-num 10 \
    --task-num 1 \
    --task-id 1 \
    --out ${section_06_dir}/vmeQTL_results/Trans_candidateSNPs/vQTL_${vQTL_method}_trans_methylation_chr${chr}
fi
echo "06d chunk${chunk} method ${vQTL_method} has been done successfully at $(date)"
