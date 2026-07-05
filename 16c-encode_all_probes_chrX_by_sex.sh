#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

mkdir -p "${section_16_dir}/logs_c"
exec &> >(tee ${section_16c_logfile})
print_version

echo "Encoding Module 16 chrX genotype and sex-specific all-probe phenotypes"
echo "This script uses the HASE mamba environment and consumes phenotypes from 16a."

mkdir -p ${hase16_chrx_encoding_female}
mkdir -p ${hase16_chrx_encoding_male}

rm -rf ${hase16_chrx_encoding_female:?}/*
rm -rf ${hase16_chrx_encoding_male:?}/*

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

    if [ ! -f "${pheno_dir}/methylation_data.csv" ]
    then
        echo "Skipping ${sex_label}: module 16 phenotype input is missing: ${pheno_dir}/methylation_data.csv"
        echo "Run 16a-prepare_chrX_all_probe_phenotypes_by_sex.sh before 16c."
        return 1
    fi

    echo "Encoding sex-specific all available methylation probes against chrX genotypes for ${sex_label} samples"
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

echo "Successfully encoded module 16 chrX genotype and sex-specific all-probe methylation data"
