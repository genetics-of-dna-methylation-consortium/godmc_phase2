#!/bin/bash

./resources/setup.sh
exec &> >(tee ${section_04e_logfile})
print_version

#Please read resources/bin/hase/README_2.md
#An example is also provided below

mkdir -p ${hase_single_site}
#mkdir -p ${hase_cov}
#awk -v OFS='\t' '{print $1,1}' <${hase_dir_in}/data.fam >${hase_cov}/covariates.txt

python ${hase}/hase.py \
   -mode single-meta \
   -study_name ${study_name} \
   -g ${hase_converting} \
   -ph ${hase_pheno} \
   -cov ${hase_cov}  \
   -mapper ${hase_mapping} \
   -o ${hase_single_site} \
   -ref_name ref-hrc

mkdir -p ${home_directory}/results/04/meta_inputs/part_dev
mkdir -p ${home_directory}/results/04/meta_inputs/use_data
mkdir -p ${home_directory}/results/04/meta_inputs/mapping
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/phenotypes
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/individuals
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/probes
mkdir -p ${home_directory}/results/04/meta_inputs/use_data/genotype

mv ${hase_single_site}/*npy ${home_directory}/results/04/meta_inputs/part_dev
mv ${hase_encoding}/encode_genotype/*h5 ${home_directory}/results/04/meta_inputs/use_data/genotype/
mv ${hase_encoding}/encode_individuals/*h5 ${home_directory}/results/04/meta_inputs/use_data/individuals/
mv ${hase_converting}/probes/$study_name.h5 ${home_directory}/results/04/meta_inputs/use_data/probes/
mv ${hase_mapping}/*npy ${home_directory}/results/04/mapping/
mv ${hase_encoding}/encode_phenotype/*.csv ${home_directory}/results/04/meta_inputs/use_data/phenotypes/

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
