#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}/logs_g
mkdir -p ${section_06_dir}/GEI_trans/candidate_CpGs

exec &> >(tee ${section_06g_logfile}_chr${1})
print_version

mamba activate tensorqtl_godmc

chr=$1

echo "Start to run 06g chr ${chr} at $(date)"
python ${scripts_directory}/resources/methylation/interaction_trans.py \
    ${tabfile}.prunedSNPs.chr${chr} \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_cpg_of_interest.bed.gz \
    ${envs_input} \
    ${section_06_dir}/GEI_trans/candidate_CpGs/snp_chr${chr} \
    0 \
    0.05

echo "06g chr ${chr} has been done successfully at $(date)"
