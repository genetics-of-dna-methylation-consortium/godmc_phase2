#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}/logs_e
mkdir -p ${section_06_dir}/vmeQTL_results/Trans_candidateCpGs/

exec &> >(tee ${section_06e_logfile}_method${1}_chr${2}_chunk${3})
print_version

vQTL_method=$1
chr=$2
chunk=$3

echo "Started to run 06e chr ${chr} chunk ${chunk} method ${vQTL_method} at $(date)"
if [ $vQTL_method = "BF" ]
then
${osca} \
    --vqtl \
    --vqtl-mtd 2 \
    --bfile ${tabfile}.prunedSNPs.chr${chr}.chunk${chunk} \
    --befile ${meth_vmeQTL_directory}/vmeQTL_phase2/cpgs_of_interest \
    --trans \
    --trans-wind 2000 \
    --thread-num 10 \
    --task-num 1 \
    --task-id 1 \
    --out ${section_06_dir}/vmeQTL_results/Trans_candidateCpGs/vQTL_${vQTL_method}_trans_SNP_chr${chr}_chunk${chunk}
else
${osca_new} \
    --vqtl \
    --vqtl-mtd ${vQTL_method} \
    --geno ${tabfile}.prunedSNPs.chr${chr}.chunk${chunk} \
    --pheno-bod ${meth_vmeQTL_directory}/vmeQTL_phase2/cpgs_of_interest \
    --trans \
    --trans-wind 2000000 \
    --thread-num 10 \
    --task-num 1 \
    --task-id 1 \
    --out ${section_06_dir}/vmeQTL_results/Trans_candidateCpGs/vQTL_${vQTL_method}_trans_SNP_chr${chr}_chunk${chunk}
fi

echo "06e chr ${chr} chunk ${chunk} method ${vQTL_method} has been done successfully at $(date)"
