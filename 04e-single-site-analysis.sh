#!/bin/bash

source resources/setup.sh "$@"
set -- $concatenated

exec &> >(tee ${section_04e_logfile})
print_version

#Please read resources/bin/hase/README_2.md
#An example is also provided below

mkdir -p ${light_hase_single_site}
#mkdir -p ${hase_cov}
#awk -v OFS='\t' '{print $1,1}' <${hase_dir_in}/data.fam >${hase_cov}/covariates.txt

"${PYTHON_RUNNER[@]}" "${light_hase}/hase.py" \
   -mode single-meta \
   -study_name ${study_name} \
   -g ${light_hase_converting} \
   -ph ${light_hase_pheno} \
   -cov ${hase_cov}  \
   -mapper ${light_hase_mapping} \
   -o ${light_hase_single_site} \
   -ref_name ref-hrc

mkdir -p ${home_directory}/results/04/meta_inputs/part_dev
mkdir -p ${home_directory}/results/04/meta_inputs/use_data
mkdir -p ${home_directory}/results/04/meta_inputs/mapping
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/phenotypes
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/individuals
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/probes
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/genotype

mv ${light_hase_single_site}/*npy ${home_directory}/results/04/meta_inputs/part_dev
mv ${light_hase_encoding}/encode_genotype/*h5 ${home_directory}/results/04/meta_inputs/use_data/genotype/
mv ${light_hase_encoding}/encode_individuals/*h5 ${home_directory}/results/04/meta_inputs/use_data/individuals/
mv ${light_hase_converting}/probes/$study_name.h5 ${home_directory}/results/04/meta_inputs/use_data/probes/
mv ${light_hase_mapping}/*npy ${home_directory}/results/04/meta_inputs/mapping/
mv ${light_hase_encoding}/encode_phenotype/*.csv ${home_directory}/results/04/meta_inputs/use_data/phenotypes/

#Example
#mkdir -p ./test/single_site_files

#${Python_directory}python hase.py \
   #-mode single-meta \
   #-study_name go_dmc_1 \
   #-g ./test/converted_files \
   #-ph ./test/hase_pheno \
   #-cov ./test/hase_cov  \
   #-mapper ./test/mapped_files  \
   #-o ./test/single_site_files \
   #-ref_name ref-hrc

echo "Single site analysis successfully completed"
