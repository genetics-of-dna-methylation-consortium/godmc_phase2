#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_05a_logfile})
print_version

# this script is to observe the female/male genotype data, and generate mapper

mkdir -p ${hase_in_female}
mkdir -p ${hase_in_male}
mkdir -p ${hase_converting_female}
mkdir -p ${hase_converting_male}
mkdir -p ${hase_mapping_female}
mkdir -p ${hase_mapping_male}

# female samples
if [ -f ${transformed_methylation_adjusted_pcs}.Female.chrX.csv ];
then
    awk -F' ' '$34 == "F" {print $1}' ${covariates_combined}.txt > ${hase_in_female}/female_id
    awk 'NR==FNR {ids[$1]; next} $2 in ids' ${hase_in_female}/female_id ${hase_dir_in}/data.fam | cut -f 1-2 > ${hase_in_female}/female_fid_id
    ${plink2} --bfile ${hase_dir_in}/data --keep ${hase_in_female}/female_fid_id --make-bed --out ${hase_in_female}/data
    rm ${hase_in_female}/data.log
    rm ${hase_in_female}/*id*

    echo "Start converting genetic data of female samples"
    python ${hase}/hase.py \
        -mode converting \
        -g ${hase_in_female} \
        -o ${hase_converting_female} \
        -study_name ${study_name} # the name for your study

    echo "Start inverting allele positions of female samples"
    python  ${hase}/added/invert_probes.py \
        -f  ${hase_converting_female}/probes \
        -n ${study_name}

    echo "Start mapping genetic data of female samples"
    python ${hase}/tools/mapper.py \
        -g ${hase_converting_female} \
        -o ${hase_mapping_female} \
        -study_name ${study_name} \
        -ref_name "ref-hrc"
else
    echo "file ${transformed_methylation_adjusted_pcs}.Female.chrX.csv does not exist, please check if no female samples in your dataset"
fi

# male samples
if [ -f ${transformed_methylation_adjusted_pcs}.Male.chrX.csv ];
then
    awk -F' ' '$34 == "M" {print $1}' ${covariates_combined}.txt > ${hase_in_male}/male_id
    awk 'NR==FNR {ids[$1]; next} $2 in ids' ${hase_in_male}/male_id ${hase_dir_in}/data.fam | cut -f 1-2 > ${hase_in_male}/male_fid_id
    ${plink2} --bfile ${hase_dir_in}/data --keep ${hase_in_male}/male_fid_id --make-bed --out ${hase_in_male}/data
    rm ${hase_in_male}/data.log
    rm ${hase_in_male}/*id*

    echo "Start converting genetic data of male samples"
    python ${hase}/hase.py \
        -mode converting \
        -g ${hase_in_male} \
        -o ${hase_converting_male} \
        -study_name ${study_name} # the name for your study 

    echo "Start inverting allele positions of male samples"
    python  ${hase}/added/invert_probes.py \
        -f  ${hase_converting_male}/probes \
        -n ${study_name}

    echo "Start mapping genetic data of male samples"
    python ${hase}/tools/mapper.py \
        -g ${hase_converting_male} \
        -o ${hase_mapping_male} \
        -study_name ${study_name} \
        -ref_name "ref-hrc"
else
    echo "file ${transformed_methylation_adjusted_pcs}.Male.chrX.csv does not exist, please check if no male samples in your dataset"
fi

echo "Successfully finished 05a script"
