#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}/logs_g
mkdir -p ${section_06_dir}/GEI_trans/candidate_CpGs

exec &> >(tee ${section_06g_logfile}_chr${1}_chunk${2})
print_version

chr=$1
chunk=$2
mkdir -p ${section_06_dir}/GEI_trans/candidate_CpGs/snp_chr${chr}/chunk${chunk}

echo "Start to run 06g chr ${chr} at $(date)"
python ${scripts_directory}/resources/methylation/interaction_trans.py \
    ${tabfile}.prunedSNPs.chr${chr}.chunk${chunk} \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_cpg_of_interest.bed.gz \
    ${envs_input} \
    ${section_06_dir}/GEI_trans/candidate_CpGs/snp_chr${chr}/chunk${chunk}/chr${chr}_chunk${chunk} \
    0 \
    0.05

${R_directory}Rscript ${scripts_directory}/resources/methylation/filter_GEI.R ${section_06_dir}/GEI_trans/candidate_CpGs/snp_chr${chr}/chunk${chunk}

if [ -f ${section_06_dir}/GEI_trans/candidate_CpGs/snp_chr${chr}/chunk${chunk}/GEI_geneticPC_interaction_5e-8.csv ];
then
    rm ${section_06_dir}/GEI_trans/candidate_CpGs/snp_chr${chr}/chunk${chunk}/*genetic_pc*parquet
fi

echo "06g chr ${chr} chunk ${chunk} has been done successfully at $(date)"
