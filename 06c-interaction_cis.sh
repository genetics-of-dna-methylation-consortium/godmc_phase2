#!/bin/bash -l

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}/logs_c
mkdir -p ${section_06_dir}/GEI_cis/

exec &> >(tee ${section_06c_logfile}_chunk${1}_chr${2})
print_version

#mamba activate tensorqtl_godmc

echo "Start to run 06c chunk${1} chr${2} at $(date)"
chunk=$1
chr=$2

mkdir -p ${section_06_dir}/GEI_cis/chunk_${chunk}_chr_${chr}

python ${scripts_directory}/resources/methylation/interaction_cis.py \
    ${vmeQTL_list1} \
    ${tabfile}.tab.${chunk} \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_cisCpGs_chr${chr}.bed.gz \
    ${chunk} \
    ${chr} \
    ${envs_input} \
    ${section_06_dir}/GEI_cis/chunk_${chunk}_chr_${chr}

${R_directory}Rscript ${scripts_directory}/resources/methylation/filter_GEI.R ${section_06_dir}/GEI_cis/chunk_${chunk}_chr_${chr}

if [ -f ${section_06_dir}/GEI_cis/chunk_${chunk}_chr_${chr}/GEI_geneticPC_interaction_5e-8.csv ];
then
    rm GEI_chunk${chunk}_chr${chr}_E_genetic_pc*candidate*parquet
fi

echo "06c chunk${1} chr${2} has been done successfully at $(date)"
