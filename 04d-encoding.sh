#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04d_logfile})
print_version

mkdir -p ${light_hase_encoding}
mkdir -p ${light_hase_pheno}
cp ${transformed_methylation_adjusted_pcs}.csv ${light_hase_pheno}
mv ${light_hase_pheno}/transformed_methylation_adjusted_pcs.csv ${light_hase_pheno}/methylation_data.csv

"${PYTHON_RUNNER[@]}" "${light_hase}/hase.py" \
   -mode encoding \
   -study_name ${study_name} \
   -g ${light_hase_converting} \
   -o ${light_hase_encoding} \
   -mapper ${light_hase_mapping} \
   -ph ${light_hase_pheno} \
   -ref_name ref-hrc

echo "Successfully encoded the genetic data"
