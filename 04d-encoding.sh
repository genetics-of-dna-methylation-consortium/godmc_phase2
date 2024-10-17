#!/bin/bash

./resources/setup.sh "$@"
exec &> >(tee ${section_04d_logfile})
print_version

#Please read resources/bin/hase/README_2.md
#An example is also provided below

mkdir -p ${hase_encoding}
mkdir -p ${hase_pheno}
cp ${transformed_methylation_adjusted_pcs}.csv ${hase_pheno}
mv ${hase_pheno}/transformed_methylation_adjusted_pcs.csv ${hase_pheno}/methylation_data.csv

python ${hase}/hase.py \
   -mode encoding \
   -study_name ${study_name} \
   -g ${hase_converting} \
   -o ${hase_encoding} \
   -mapper ${hase_mapping} \
   -ph ${hase_pheno} \
   -ref_name ref-hrc

echo "Successfully encoded the genetic data"
