#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_05b_logfile})
print_version

#Please read resources/bin/light_hase/README_2.md
#An example is also provided below

mkdir -p ${hase_encoding_female}
mkdir -p ${hase_pheno_female}
mkdir -p ${hase_encoding_male}
mkdir -p ${hase_pheno_male}

rm -rf ${hase_encoding_female:?}/*
rm -rf ${hase_pheno_female:?}/*
rm -rf ${hase_encoding_male:?}/*
rm -rf ${hase_pheno_male:?}/*

# add if file exist check - in case all female / all male
if [ -f ${transformed_methylation_adjusted_pcs}.Female.chrX.csv ];
then
    cp ${transformed_methylation_adjusted_pcs}.Female.chrX.csv ${hase_pheno_female}
    mv ${hase_pheno_female}/transformed_methylation_adjusted_pcs.Female.chrX.csv ${hase_pheno_female}/methylation_data.csv
    python ${light_hase}/hase.py \
        -mode encoding \
        -study_name ${study_name} \
        -g ${hase_converting_female} \
        -o ${hase_encoding_female} \
        -mapper ${hase_mapping_female} \
        -ph ${hase_pheno_female} \
        -ref_name ref-hrc
else
    echo "file ${transformed_methylation_adjusted_pcs}.Female.chrX.csv does not exist, please check if no female samples in your dataset"
fi

if [ -f ${transformed_methylation_adjusted_pcs}.Male.chrX.csv ];
then
    if [ ! -f ${transformed_methylation_adjusted_pcs}.Male.chrY.csv ];
    then
        echo "ERROR: Missing male chrY methylation file: ${transformed_methylation_adjusted_pcs}.Male.chrY.csv"
        echo "Male 05b encoding requires both male chrX and chrY files to create the combined phenotype input"
        exit 1
    fi
    cat ${transformed_methylation_adjusted_pcs}.Male.chrX.csv <(tail -n +2 ${transformed_methylation_adjusted_pcs}.Male.chrY.csv) > ${transformed_methylation_adjusted_pcs}.Male.chrX.chrY.csv
    cp ${transformed_methylation_adjusted_pcs}.Male.chrX.chrY.csv ${hase_pheno_male}
    mv ${hase_pheno_male}/transformed_methylation_adjusted_pcs.Male.chrX.chrY.csv ${hase_pheno_male}/methylation_data.csv
    python ${light_hase}/hase.py \
        -mode encoding \
        -study_name ${study_name} \
        -g ${hase_converting_male} \
        -o ${hase_encoding_male} \
        -mapper ${hase_mapping_male} \
        -ph ${hase_pheno_male} \
        -ref_name ref-hrc
else
    echo "file ${transformed_methylation_adjusted_pcs}.Male.chrX.csv does not exist, please check if no male samples in your dataset"
fi

echo "Successfully encoded the genetic data"
