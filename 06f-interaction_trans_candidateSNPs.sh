#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_06_dir}/logs_f
mkdir -p ${section_06_dir}/GEI_trans/candidate_SNPs/cpg_chunk${chunk}
mkdir -p ${section_06_dir}/GEI_trans/epistasis

exec &> >(tee ${section_06f_logfile}_chunk${1})
print_version

chunk=$1

echo "Start to run 06f CpG chunk ${chunk} at $(date)"
python ${scripts_directory}/resources/methylation/interaction_trans.py \
    ${tabfile}.vQTLs \
    ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_allCpGs_chunk${chunk}.bed.gz \
    ${envs_input} \
    ${section_06_dir}/GEI_trans/candidate_SNPs/cpg_chunk${chunk}/GEI_trans_cpg_chunk${chunk} \
    0 \
    0.05

for row in $(seq 1 5);
do
    if [ -f ${tabfile}_epi_row${row}.raw ];
    then
        awk 'BEGIN{FS=" ";OFS="\t"}{print $2,$7}' ${tabfile}_epi_row${row}.raw > ${tabfile}_epi_row${row}.raw1
        grep -v NA ${tabfile}_epi_row${row}.raw1 > ${tabfile}_epi_row${row}.raw2
        python ${scripts_directory}/resources/methylation/interaction_trans.py \
            ${tabfile}_epi_row${row} \
            ${meth_vmeQTL_directory}/vmeQTL_phase2/adjustcovs_cpg_phase2_allCpGs_chunk${chunk}.bed.gz \
            ${tabfile}_epi_row${row}.raw2 \
            ${section_06_dir}/GEI_trans/epistasis/epistasis_cpg_chunk${chunk} \
            0 \
            1
    else
        echo "Skipping epistasis pair $row"
    fi
done
echo "06f CpG chunk ${chunk} has been done at $(date)"
