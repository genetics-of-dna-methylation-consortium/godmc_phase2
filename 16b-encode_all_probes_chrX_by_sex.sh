#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_16b_logfile})
print_version

mkdir -p ${hase16_chrx_encoding_female}
mkdir -p ${hase16_chrx_encoding_male}
mkdir -p ${hase16_allprobes_pheno_female}
mkdir -p ${hase16_allprobes_pheno_male}

rm -rf ${hase16_chrx_encoding_female:?}/*
rm -rf ${hase16_chrx_encoding_male:?}/*
rm -rf ${hase16_allprobes_pheno_female:?}/*
rm -rf ${hase16_allprobes_pheno_male:?}/*

if [ ! -f "${transformed_methylation_adjusted_pcs}.csv" ]
then
    echo "ERROR: Missing all-probe methylation file: ${transformed_methylation_adjusted_pcs}.csv"
    exit 1
fi

encode_chrx_all_probes() {
    sex_label="$1"
    converting_dir="$2"
    mapping_dir="$3"
    pheno_dir="$4"
    encoding_dir="$5"

    if [ ! -d "${converting_dir}/probes" ]
    then
        echo "Skipping ${sex_label}: chrX converted genotype directory is missing: ${converting_dir}/probes"
        return 1
    fi

    if [ ! -d "${mapping_dir}" ]
    then
        echo "Skipping ${sex_label}: chrX mapping directory is missing: ${mapping_dir}"
        return 1
    fi

    cp ${transformed_methylation_adjusted_pcs}.csv ${pheno_dir}
    mv ${pheno_dir}/transformed_methylation_adjusted_pcs.csv ${pheno_dir}/methylation_data.csv

    echo "Encoding all methylation probes against chrX genotypes for ${sex_label} samples"
    python ${light_hase}/hase.py \
        -mode encoding \
        -study_name ${study_name} \
        -g ${converting_dir} \
        -o ${encoding_dir} \
        -mapper ${mapping_dir} \
        -ph ${pheno_dir} \
        -ref_name ref-hrc
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase encoding failed for ${sex_label} module 16 inputs"
        exit 1
    fi

    return 0
}

encoded_sex_count=0

if encode_chrx_all_probes \
    "female" \
    "${hase16_chrx_converting_female}" \
    "${hase16_chrx_mapping_female}" \
    "${hase16_allprobes_pheno_female}" \
    "${hase16_chrx_encoding_female}"
then
    encoded_sex_count=$((encoded_sex_count + 1))
fi

if encode_chrx_all_probes \
    "male" \
    "${hase16_chrx_converting_male}" \
    "${hase16_chrx_mapping_male}" \
    "${hase16_allprobes_pheno_male}" \
    "${hase16_chrx_encoding_male}"
then
    encoded_sex_count=$((encoded_sex_count + 1))
fi

if [ "${encoded_sex_count}" -eq "0" ]
then
    echo "ERROR: No female or male module 16 inputs were encoded"
    exit 1
fi

echo "Successfully encoded module 16 chrX genotype and all-probe methylation data"
