#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_16_dir}/logs_a"
exec &> >(tee ${section_16a_logfile})
print_version

mkdir -p ${hase16_allprobes_pheno_female}
mkdir -p ${hase16_allprobes_pheno_male}

rm -f "${hase16_allprobes_pheno_female}/methylation_data.csv"
rm -f "${hase16_allprobes_pheno_male}/methylation_data.csv"

echo "Preparing Module 16 sex-specific all-probe methylation phenotypes"
echo "This script uses the R environment. Activate the HASE mamba environment after 16a."

${R_directory}Rscript resources/methylation/prepare_module16_phenotypes.R \
    "${covariates_combined}.txt" \
    "${bfile}.fam" \
    "${transformed_methylation_adjusted_pcs}" \
    "${hase16_allprobes_pheno_female}" \
    "${hase16_allprobes_pheno_male}"

if [ "$?" -ne "0" ]
then
    echo "ERROR: Module 16 phenotype preparation failed"
    exit 1
fi

echo "Module 16 sex-specific all-probe methylation phenotypes successfully prepared"
