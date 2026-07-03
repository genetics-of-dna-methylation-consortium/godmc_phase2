#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_16c_logfile})
print_version

mkdir -p ${hase16_chrx_single_site_female}
mkdir -p ${hase16_chrx_single_site_male}

rm -rf ${hase16_chrx_single_site_female:?}/*
rm -rf ${hase16_chrx_single_site_male:?}/*

prepare_meta_input_dirs() {
    meta_inputs="$1"

    mkdir -p ${meta_inputs}/part_dev
    mkdir -p ${meta_inputs}/use_data
    mkdir -p ${meta_inputs}/mapping
    mkdir -p ${meta_inputs}/use_data/phenotypes
    mkdir -p ${meta_inputs}/use_data/individuals
    mkdir -p ${meta_inputs}/use_data/probes
    mkdir -p ${meta_inputs}/use_data/genotype

    rm -rf ${meta_inputs}/part_dev/*
    rm -rf ${meta_inputs}/mapping/*
    rm -rf ${meta_inputs}/use_data/phenotypes/*
    rm -rf ${meta_inputs}/use_data/individuals/*
    rm -rf ${meta_inputs}/use_data/probes/*
    rm -rf ${meta_inputs}/use_data/genotype/*
}

run_chrx_all_probe_single_meta() {
    sex_label="$1"
    converting_dir="$2"
    mapping_dir="$3"
    pheno_dir="$4"
    encoding_dir="$5"
    single_site_dir="$6"
    meta_inputs="$7"

    if [ ! -f "${pheno_dir}/methylation_data.csv" ]
    then
        echo "Skipping ${sex_label}: module 16 phenotype input is missing: ${pheno_dir}/methylation_data.csv"
        return 1
    fi

    if [ ! -d "${encoding_dir}/encode_genotype" ]
    then
        echo "Skipping ${sex_label}: encoded genotype directory is missing: ${encoding_dir}/encode_genotype"
        return 1
    fi

    echo "Running module 16 chrX sex-specific all-probe HASE single-meta for ${sex_label} samples"
    python ${light_hase}/hase.py \
        -mode single-meta \
        -study_name ${study_name} \
        -g ${converting_dir} \
        -ph ${pheno_dir} \
        -cov ${hase_cov}  \
        -mapper ${mapping_dir} \
        -o ${single_site_dir} \
        -ref_name ref-hrc
    if [ "$?" -ne "0" ]
    then
        echo "ERROR: light_hase single-meta failed for ${sex_label} module 16 inputs"
        exit 1
    fi

    prepare_meta_input_dirs "${meta_inputs}"

    cp ${single_site_dir}/*npy ${meta_inputs}/part_dev
    cp ${encoding_dir}/encode_genotype/*h5 ${meta_inputs}/use_data/genotype/
    cp ${encoding_dir}/encode_individuals/*h5 ${meta_inputs}/use_data/individuals/
    cp ${converting_dir}/probes/$study_name.h5 ${meta_inputs}/use_data/probes/
    cp ${mapping_dir}/*npy ${meta_inputs}/mapping/
    cp ${encoding_dir}/encode_phenotype/*.csv ${meta_inputs}/use_data/phenotypes/

    return 0
}

single_meta_sex_count=0

if run_chrx_all_probe_single_meta \
    "female" \
    "${hase16_chrx_converting_female}" \
    "${hase16_chrx_mapping_female}" \
    "${hase16_allprobes_pheno_female}" \
    "${hase16_chrx_encoding_female}" \
    "${hase16_chrx_single_site_female}" \
    "${section_16_dir}/meta_inputs_female"
then
    single_meta_sex_count=$((single_meta_sex_count + 1))
fi

if run_chrx_all_probe_single_meta \
    "male" \
    "${hase16_chrx_converting_male}" \
    "${hase16_chrx_mapping_male}" \
    "${hase16_allprobes_pheno_male}" \
    "${hase16_chrx_encoding_male}" \
    "${hase16_chrx_single_site_male}" \
    "${section_16_dir}/meta_inputs_male"
then
    single_meta_sex_count=$((single_meta_sex_count + 1))
fi

if [ "${single_meta_sex_count}" -eq "0" ]
then
    echo "ERROR: No female or male module 16 single-meta analyses were completed"
    exit 1
fi

echo "Module 16 chrX genotype against sex-specific all methylation probes successfully completed"
