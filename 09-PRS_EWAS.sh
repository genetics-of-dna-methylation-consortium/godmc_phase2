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

####
# additional glint section for testing

if["${related}" = "yes"]
then

# get unrelated IDs
echo "Removing any cryptic relateds"

${gcta} \
    --grm ${grmfile_all} \
    --grm-cutoff ${rel_cutoff} \
    --make-grm-bin \
    --out ${grmfile_glint_unrelated}

mamba activate deep_env

# convert grmfile_all to txt file - all participants 
echo "convert grm to txt file"
  ${R_directory}Rscript resources/ewas/format_data.R \
    ${grm_base} \
    ${methylation_no_outliers} \
    ADHD \ # changed from ${PRS} because we only need ADHD
    ${nongenetic_meth_pcs_untransformed} \
    ${DEEP_scripts_directory} \
    ${glint_output_path} \
    all

  echo "Successfully completed conversion to text files for glint (all participants)"

# convert grmfile_all to txt file - unrelated participants 
echo "convert grm to txt file"
  ${R_directory}Rscript resources/ewas/format_data.R \
    ${grm_glint_unrelated_base} \
    ${methylation_no_outliers} \
    ADHD \ # changed from ${PRS} because we only need ADHD
    ${nongenetic_meth_pcs_untransformed} \
    ${DEEP_scripts_directory} \
    ${glint_output_path} \
    unrelated

  echo "Successfully completed conversion to text files for glint (unrelated participants)"


mamba activate hase_py2

echo "convert data to glint format"
 
# convert data to glint format
python ${glint_directory}/glint-1.0.4/glint.py --datafile  ${glint_output_path}/dnam_for_glint_all.txt --covarfile /${glint_output_path}/covariates_for_glint_all.txt --phenofile /${glint_output_path}/phenotypes_for_glint_all.txt --gsave --out ${glint_output_path}/datafile_for_glint_all

python ${glint_directory}/glint-1.0.4/glint.py --datafile  ${glint_output_path}/dnam_for_glint_unrelated.txt --covarfile /${glint_output_path}/covariates_for_glint_unrelated.txt --phenofile /${glint_output_path}/phenotypes_for_glint_unrelated.txt --gsave --out ${glint_output_path}/datafile_for_glint_unrelated


  echo "Successfully completed conversion to glint format"

echo "run glint ewas"

# run glint EWAS
python ${glint_directory}/glint-1.0.4/glint.py --datafile ${glint_output_path}/datafile_for_glint_all.glint --ewas --lmm --pheno ADHD --kinship ${glint_output_path}/grm_for_glint_all.txt  --out ${glint_output_path}/glint_ewas_all

python ${glint_directory}/glint-1.0.4/glint.py --datafile ${glint_output_path}/datafile_for_glint_unrelated.glint --ewas --lmm --pheno ADHD --kinship ${glint_output_path}/grm_for_glint_unrelated.txt  --out ${glint_output_path}/glint_ewas_unrelated


  echo "Successfully completed glint EWAS"

mamba activate deep_env

# run plots script
echo "running glint plots"
  ${R_directory}Rscript resources/ewas/glint_plots.R \
    ${glint_ewas_all} \
    ${glint_output_path} \
    ${study_name} \
    ${meth_array} \
    all

  ${R_directory}Rscript resources/ewas/glint_plots.R \
    ${glint_ewas_unrelated} \
    ${glint_output_path} \
    ${study_name} \
    ${meth_array} \
    unrelated


  echo "Successfully completed glint Manhattan and qq plots"

# run glint-original ewas comparison
echo "running glint comparison"
  ${R_directory}Rscript resources/ewas/original_glint_comparison.R \
    ${glint_ewas_all} \
    ${glint_output_path} \
    ${section_09_dir}/ADHD \
    ${study_name} \
    ADHD \
    all

  ${R_directory}Rscript resources/ewas/original_glint_comparison.R \
    ${glint_ewas_unrelated} \
    ${glint_output_path} \
    ${section_09_dir}/ADHD \
    ${study_name} \
    ADHD \
    unrelated


  echo "Successfully completed glint-original ewas comparison"

# run all-unrelated glint ewas comparison
echo "running all-unrelated glint comparison"
  ${R_directory}Rscript resources/ewas/glint_plots_all_unrelated.R \
    ${glint_ewas_all} \
    ${glint_ewas_unrelated} \
    ${glint_output_path} \
    ${section_09_dir}/ADHD \
    ${study_name} \
    ADHD

  echo "Successfully completed all-unrelated ewas comparison"


# end of additional glint section
fi

####

 
} | tee "$log_file"
  
done

echo ""
echo "script finalised"
echo ""
 
