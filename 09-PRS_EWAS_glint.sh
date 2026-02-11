#!/bin/bash

set -e -o pipefail

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p ${section_09_dir}/logs
exec &> >(tee ${section_09_logfile})
print_version

####
# additional glint section for testing
#source ~/.bashrc
#conda activate hase_py2_glint

if [ "${related}" = "yes" ]
then
mkdir -p $glint_output_path
# get unrelated IDs
echo "Removing any cryptic relateds"

${gcta} \
    --grm ${grmfile_all} \
    --grm-cutoff ${rel_cutoff} \
    --make-grm-bin \
    --out ${grmfile_glint_unrelated}

echo "convert grm to txt file"
# changed from ${PRS} because we only need ADHD
  ${R_directory}Rscript resources/ewas/format_data.R \
    "${grm_base}" \
    "${methylation_no_outliers}" \
    "ADHD" \
    "${covariates_combined}.txt" \
    "${glint_output_path}" \
    "all" \
    "${home_directory}"

  echo "Successfully completed conversion to text files for glint (all participants)"

# convert grmfile_all to txt file - unrelated participants 
# changed from ${PRS} because we only need ADHD
echo "convert grm to txt file"
  ${R_directory}Rscript resources/ewas/format_data.R \
    "${grm_glint_unrelated_base}" \
    "${methylation_no_outliers}" \
    "ADHD" \
    "${covariates_combined}.txt" \
    "${glint_output_path}" \
    "unrelated" \
    "${home_directory}"

echo "Successfully completed conversion to text files for glint (unrelated participants)"

echo "convert data to glint format"
 
# convert data to glint format
$glint --datafile  ${glint_output_path}/dnam_for_glint_all.txt --covarfile /${glint_output_path}/covariates_for_glint_all.txt --phenofile /${glint_output_path}/phenotypes_for_glint_all.txt --gsave --out ${glint_output_path}/datafile_for_glint_all

$glint --datafile  ${glint_output_path}/dnam_for_glint_unrelated.txt --covarfile /${glint_output_path}/covariates_for_glint_unrelated.txt --phenofile /${glint_output_path}/phenotypes_for_glint_unrelated.txt --gsave --out ${glint_output_path}/datafile_for_glint_unrelated


  echo "Successfully completed conversion to glint format"

echo "run glint ewas"

# run glint EWAS
$glint --datafile ${glint_output_path}/datafile_for_glint_all.glint --ewas --lmm --pheno ADHD --kinship ${glint_output_path}/grm_for_glint_all.txt  --out ${glint_output_path}/glint_ewas_all

$glint --datafile ${glint_output_path}/datafile_for_glint_unrelated.glint --ewas --lmm --pheno ADHD --kinship ${glint_output_path}/grm_for_glint_unrelated.txt  --out ${glint_output_path}/glint_ewas_unrelated


  echo "Successfully completed glint EWAS"

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
 
