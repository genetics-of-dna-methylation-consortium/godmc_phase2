#!/bin/bash

set -e -o pipefail

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_09_dir}/logs
exec &> >(tee ${section_09_logfile})
print_version

vect_PRS=$(grep "PRS" ${scripts_directory}/resources/parameters | grep "weights" | awk -F"_" '{print $2}' |tr "\n" " ")
vect_PRS_array=($vect_PRS)

vect_PRS_weights=$(grep "PRS" ${scripts_directory}/resources/parameters | grep "weights" | awk -F"=" '{print $1}'  |tr "\n" " ")
vect_PRS_weights_array=($vect_PRS_weights)

n=$((${#vect_PRS_weights_array[*]}-1))

for ((k=0;k<=$n;k++))
do

PRS=${vect_PRS_array[$k]}

log_dir=${section_09_dir}/${PRS}/logs
log_file=${section_09_dir}/${PRS}/logs/log.txt

mkdir -p $log_dir

{
print_version

PRS_file=${home_directory}/processed_data/genetic_data/PRS_${PRS}
PRS_weights=${vect_PRS_weights_array[$k]}
pheno_for_PRS=phenotypes_${PRS}
cov_for_PRS=covariates_${PRS}

echo ""
echo "Generating PRS for ${PRS}"
echo ""

${plink2} \
  --bfile ${bfile} \
  --score ${!PRS_weights} 2 4 6 'list-variants' \
  --out ${PRS_file}

echo ""
echo "Standarising PRS and generating QC plots for ${PRS}"
echo ""

${R_directory}Rscript ${scripts_directory}/resources/genetics/PRS_qc.R \
  ${PRS} \
  ${PRS_file}.sscore \
  ${section_09_dir}/${PRS} \
  ${!pheno_for_PRS} \
  ${cellcounts_cov} \
  ${nongenetic_meth_pcs_untransformed} \
  ${study_name}

echo ""
echo "Running EWAS for ${PRS}"
echo ""


${R_directory}Rscript ${scripts_directory}/resources/methylation/PRS.ewas.meffil.R \
  ${untransformed_methylation_adjusted}.RData \
  ${PRS} \
  ${!cov_for_PRS} \
  ${nongenetic_meth_pcs_untransformed} \
  ${home_directory} \
  ${section_09_dir}/${PRS} \
  ${study_name}

echo ""
echo "EWAS run successfully for ${PRS}"
echo ""
 
} | tee "$log_file"
  
done

echo ""
echo "script finalised"
echo ""
 
