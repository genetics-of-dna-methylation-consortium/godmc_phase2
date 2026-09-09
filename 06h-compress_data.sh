#!/bin/bash
source resources/setup.sh "$@"
set -- $concatenated

suff="tgz"
flags="czf"

#mkdir -p ${section_06_dir}/upload/results_part1
#mkdir -p ${section_06_dir}/upload/results_part2
#mkdir -p ${section_06_dir}/upload/results_part3
#mkdir -p ${section_06_dir}/upload/results_part4
#mkdir -p ${section_06_dir}/upload/results_part5
#mkdir -p ${section_06_dir}/upload/results_part6
#mkdir -p ${section_06_dir}/upload/results_part7

#mv ${section_06_dir}/GEI_cis ${section_06_dir}/upload/results_part1
#mv ${section_06_dir}/vmeQTL_results/Trans_candidateCpGs ${section_06_dir}/upload/results_part2
#mv ${section_06_dir}/vmeQTL_results/Trans_candidateSNPs ${section_06_dir}/upload/results_part3
#mv ${section_06_dir}/GEI_trans/candidate_CpGs ${section_06_dir}/upload/results_part4
#mv ${section_06_dir}/GEI_trans/candidate_SNPs ${section_06_dir}/upload/results_part5
#mv ${section_06_dir}/GEI_trans/epistasis ${section_06_dir}/upload/results_part6
#mv ${section_06_dir}/logs* ${section_06_dir}/upload/results_part7
#mv ${section_06_dir}/E_plots.pdf ${section_06_dir}/upload/results_part7
#mv ${section_06_dir}/E_summary.csv ${section_06_dir}/upload/results_part7

tar ${flags} ${home_directory}/results/${study_name}_06_results1.${suff} -C ${home_directory} ${section_06_dir}/GEI_cis
tar ${flags} ${home_directory}/results/${study_name}_06_results2.${suff} -C ${home_directory} ${section_06_dir}/vmeQTL_results/Trans_candidateCpGs
tar ${flags} ${home_directory}/results/${study_name}_06_results3.${suff} -C ${home_directory} ${section_06_dir}/vmeQTL_results/Trans_candidateSNPs
tar ${flags} ${home_directory}/results/${study_name}_06_results4.${suff} -C ${home_directory} ${section_06_dir}/GEI_trans/candidate_CpGs
tar ${flags} ${home_directory}/results/${study_name}_06_results5.${suff} -C ${home_directory} ${section_06_dir}/GEI_trans/candidate_SNPs

if [[ -d "${section_06_dir}/vmeQTL_results/Missing_association" ]]; then
    tar ${flags} ${home_directory}/results/${study_name}_06_results6.${suff} -C ${home_directory} ${section_06_dir}/GEI_trans/epistasis ${section_06_dir}/vmeQTL_results/Missing_association
else
    tar ${flags} ${home_directory}/results/${study_name}_06_results6.${suff} -C ${home_directory} ${section_06_dir}/GEI_trans/epistasis
fi

tar ${flags} ${home_directory}/results/${study_name}_06_results7.${suff} -C ${home_directory} ${section_06_dir}/logs_a ${section_06_dir}/logs_b ${section_06_dir}/logs_c ${section_06_dir}/logs_d ${section_06_dir}/logs_e ${section_06_dir}/logs_f ${section_06_dir}/logs_g ${section_06_dir}/E_plots.pdf ${section_06_dir}/E_summary.csv
