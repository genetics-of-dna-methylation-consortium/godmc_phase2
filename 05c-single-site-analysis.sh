#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_05c_logfile})
print_version

#Please read resources/bin/hase/README_2.md
#An example is also provided below

mkdir -p ${hase_single_site_female}
mkdir -p ${hase_single_site_male}

if [ -f ${transformed_methylation_adjusted_pcs}.Female.chrX.csv ];
then
    python ${hase}/hase.py \
        -mode single-meta \
        -study_name ${study_name} \
        -g ${hase_converting_female} \
        -ph ${hase_pheno_female} \
        -cov ${hase_cov}  \
        -mapper ${hase_mapping_female} \
        -o ${hase_single_site_female} \
        -ref_name ref-hrc

    mkdir -p ${home_directory}/results/05/meta_inputs_female/part_dev
    mkdir -p ${home_directory}/results/05/meta_inputs_female/use_data
    mkdir -p ${home_directory}/results/05/meta_inputs_female/mapping
    mkdir -p ${home_directory}/results/05/meta_inputs_female/use_data/phenotypes
    mkdir -p ${home_directory}/results/05/meta_inputs_female/use_data/individuals
    mkdir -p ${home_directory}/results/05/meta_inputs_female/use_data/probes
    mkdir -p ${home_directory}/results/05/meta_inputs_female/use_data/genotype

    mv ${hase_single_site_female}/*npy ${home_directory}/results/05/meta_inputs_female/part_dev
    mv ${hase_encoding_female}/encode_genotype/*h5 ${home_directory}/results/05/meta_inputs_female/use_data/genotype/
    mv ${hase_encoding_female}/encode_individuals/*h5 ${home_directory}/results/05/meta_inputs_female/use_data/individuals/
    mv ${hase_converting_female}/probes/$study_name.h5 ${home_directory}/results/05/meta_inputs_female/use_data/probes/
    mv ${hase_mapping_female}/*npy ${home_directory}/results/05/meta_inputs_female/mapping/
    mv ${hase_encoding_female}/encode_phenotype/*.csv ${home_directory}/results/05/meta_inputs_female/use_data/phenotypes/
fi

if [ -f ${transformed_methylation_adjusted_pcs}.Male.chrX.csv ];
then
    python ${hase}/hase.py \
        -mode single-meta \
        -study_name ${study_name} \
        -g ${hase_converting_male} \
        -ph ${hase_pheno_male} \
        -cov ${hase_cov}  \
        -mapper ${hase_mapping_male} \
        -o ${hase_single_site_male} \
        -ref_name ref-hrc

    mkdir -p ${home_directory}/results/05/meta_inputs_male/part_dev
    mkdir -p ${home_directory}/results/05/meta_inputs_male/use_data
    mkdir -p ${home_directory}/results/05/meta_inputs_male/mapping
    mkdir -p ${home_directory}/results/05/meta_inputs_male/use_data/phenotypes
    mkdir -p ${home_directory}/results/05/meta_inputs_male/use_data/individuals
    mkdir -p ${home_directory}/results/05/meta_inputs_male/use_data/probes
    mkdir -p ${home_directory}/results/05/meta_inputs_male/use_data/genotype

    mv ${hase_single_site_female}/*npy ${home_directory}/results/05/meta_inputs_male/part_dev
    mv ${hase_encoding_female}/encode_genotype/*h5 ${home_directory}/results/05/meta_inputs_male/use_data/genotype/
    mv ${hase_encoding_female}/encode_individuals/*h5 ${home_directory}/results/05/meta_inputs_male/use_data/individuals/
    mv ${hase_converting_female}/probes/$study_name.h5 ${home_directory}/results/05/meta_inputs_male/use_data/probes/
    mv ${hase_mapping_female}/*npy ${home_directory}/results/05/meta_inputs_male/mapping/
    mv ${hase_encoding_female}/encode_phenotype/*.csv ${home_directory}/results/05/meta_inputs_male/use_data/phenotypes/
fi

echo "Single site analysis successfully completed"
