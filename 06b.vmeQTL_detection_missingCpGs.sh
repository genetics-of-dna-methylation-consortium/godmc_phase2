#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}/logs_b
mkdir -p ${section_06_dir}/vmeQTL_results/Missing_association

exec &> >(tee ${section_06b_logfile}_chr${1})
print_version

c=$1

participate07=`grep ${study_name} ${scripts_directory}/resources/methylation/vmeQTL/vmeQTL_phase1_cohort_list.txt | wc -l`
if [ ${participate07} -eq 1 ];
then    
    echo "Started to run 06b chr${c} at $(date)"

    ${osca_new} \
    --vqtl \
    --vqtl-mtd drm \
    --geno ${tabfile}.chr${c} \
    --pheno-bod ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c} \
    --cis \
    --cis-wind 2000000 \
    --thread-num 10 \
    --task-num 1 \
    --task-id 1 \
    --out ${section_06_dir}/vmeQTL_results/Missing_association/vQTL_drm_chr${c}_missing

    ${osca_new} \
    --vqtl \
    --vqtl-mtd svlm \
    --geno ${tabfile}.chr${c} \
    --pheno-bod ${meth_vmeQTL_directory}/vmeQTL_phase2/missing_cpgs_chr${c} \
    --cis \
    --cis-wind 2000000 \
    --thread-num 10 \
    --task-num 1 \
    --task-id 1 \
    --out ${section_06_dir}/vmeQTL_results/Missing_association/vQTL_svlm_chr${c}_missing

    echo "06b chr${c} has been done successfully at $(date)"
else
    echo "06b is a supplmentary script for module 07. Your cohort was not involved in module 07 meta-analysis, please skip 06b. If you have any doubts, please contact Xiaopu (xiaopu.1.zhang@kcl.ac.uk)"
    exit
fi
